import logging
import os
import re
import time
import sys
import itertools

import click
import click_log
from tqdm import tqdm
import pickle
import enum
import operator

from polyleven import levenshtein

import pysam

from construct import *

from collections import defaultdict
from collections import Counter

import longbow.utils.constants
from ..utils import bam_utils
from ..utils.bam_utils import SegmentInfo
from ..utils.cli_utils import get_field_count_and_percent_string

from ..utils.constants import FFORMAT

PROG_NAME = "correct_umi"

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(PROG_NAME)
click_log.basic_config(logger)

MAS_GENE_PREFIX = "MAS"
GENCODE_GENE_PREFIX = "ENSG"


class ReadType(enum.Enum):
    CCS = enum.auto()
    CLR = enum.auto()


CODING_REGION = "cDNA"  # TODO: Make this take in / read in a model so that we don't have to specify this.
READ_QUALITY_TAG = "rq"


class ReadSnapshot:

    def __init__(self, read, pre_extracted, umi_tag) -> None:
        self.umi = read.get_tag(umi_tag)
        self.type = get_read_type(read)
        self.start = read.reference_start
        self.end = read.reference_end
        sequence = get_read_seq(read, pre_extracted)
        self.len = len(sequence)
        self.gc = float(sequence.count('C') + sequence.count('G'))/len(sequence)
        self.name = read.qname

    def __eq__(self, other) -> bool:
        return self.umi == other.umi and self.type == other.type and self.start == other.start and \
               self.end == other.end and self.len == other.len and \
               self.gc == other.gc and self.name == other .name

    def __hash__(self) -> int:
        return super().__hash__()


############################################################


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "--cache-read-loci",
    is_flag=True,
    default=False,
    show_default=True,
    help="Dump and read from a cache for the read loci, if possible."
)
@click.option(
    "--max-ccs-edit-dist",
    type=int,
    default=2,
    show_default=True,
    help="Maximum edit distance between a CCS read and a target to consider them both part of a group of UMIs",
)
@click.option(
    "--max-clr-edit-dist",
    type=int,
    default=3,
    show_default=True,
    help="Maximum edit distance between a CLR read and a target to consider them both part of a group of UMIs",
)
@click.option(
    "--max-ccs-length-diff",
    type=int,
    default=50,
    show_default=True,
    help="Maximum length difference between a CCS read and a target to consider them both part of a group of UMIs",
)
@click.option(
    "--max-clr-length-diff",
    type=int,
    default=150,
    show_default=True,
    help="Maximum length difference between a CLR read and a target to consider them both part of a group of UMIs",
)
@click.option(
    "--max-ccs-gc-diff",
    type=float,
    default=0.05,
    show_default=True,
    help="Maximum GC content difference between a CCS read and a target to consider them both part of a group of UMIs",
)
@click.option(
    "--max-clr-gc-diff",
    type=float,
    default=0.15,
    show_default=True,
    help="Maximum GC content difference between a CLR read and a target to consider them both part of a group of UMIs",
)
@click.option(
    "--max-ccs-umi-length-delta",
    type=int,
    default=3,
    show_default=True,
    help="Maximum length difference between the UMI of a CCS read and the given UMI length to be included in the loci "
         "for possible correction.",
)
@click.option(
    "--max-clr-umi-length-delta",
    type=int,
    default=4,
    show_default=True,
    help="Maximum length difference between the UMI of a CLR read and the given UMI length to be included in the loci "
         "for possible correction.",
)
@click.option(
    "--max-final-ccs-umi-length-delta",
    type=int,
    default=3,
    show_default=True,
    help="Maximum length difference between the final, corrected UMI of a CCS read and the given UMI length to be "
         "included in the final good results file.",
)
@click.option(
    "--max-final-clr-umi-length-delta",
    type=int,
    default=3,
    show_default=True,
    help="Maximum length difference between the final, corrected UMI of a CLR read and the given UMI length to be "
         "included in the final good results file.",
)
@click.option(
    "--min-back-seg-score",
    type=int,
    default=10,
    show_default=True,
    help="Minimum score of the back alignment tag for a read to be included in the final good results file."
)
@click.option(
    "--umi-tag",
    type=str,
    default="JX",
    show_default=True,
    help="Tag from which to read in UMIs from the input bam file."
)
@click.option(
    "--gene-tag",
    type=str,
    default="XG",
    show_default=True,
    help="Tag from which to read in gene IDs from the input bam file."
)
@click.option(
    "--eq-class-tag",
    type=str,
    default="eq",
    show_default=True,
    help="Tag from which to read in read transcript equivalence classes (i.e. transcript IDs) from the input bam file."
)
@click.option(
    "--final-umi-tag",
    type=str,
    default="BX",
    show_default=True,
    help="Tag into which to put final, corrected UMI values."
)
@click.option(
    "--umi-corrected-tag",
    type=str,
    default="UX",
    show_default=True,
    help="Tag into which to put whether a given UMI was actually corrected."
)
@click.option(
    "--back-alignment-score-tag",
    type=str,
    default="JB",
    show_default=True,
    help="Tag containing the back (trailing adapter) alignment score.",
)
@click.option(
    "-l",
    "--umi-length",
    type=int,
    default=10,
    show_default=True,
    help="Length of the UMI for this sample.",
)
@click.option(
    "-o",
    "--output-bam",
    default="-",
    type=click.Path(exists=False),
    help="Corrected UMI bam output [default: stdout].",
)
@click.option(
    "-x",
    "--reject-bam",
    type=click.Path(exists=False),
    help="Filtered bam output (failing reads only).",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.option('--pre-extracted',
              is_flag=True,
              default=False,
              show_default=True,
              help='Whether the input file has been processed with `longbow extract`')
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(umi_length, max_ccs_edit_dist, max_clr_edit_dist, max_ccs_length_diff, max_clr_length_diff, max_ccs_gc_diff,
         max_clr_gc_diff, max_ccs_umi_length_delta, max_clr_umi_length_delta, max_final_ccs_umi_length_delta,
         max_final_clr_umi_length_delta, min_back_seg_score, umi_tag, gene_tag, eq_class_tag, final_umi_tag,
         umi_corrected_tag, back_alignment_score_tag, output_bam, reject_bam, force, pre_extracted, cache_read_loci,
         input_bam):
    """Correct UMIs with Set Cover algorithm."""
    # This algorithm was originally developed by Victoria Popic and imported into Longbow by Jonn Smith.

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files([output_bam, reject_bam], exist_ok=force)

    # Load number of reads, if pbi exists:
    pbi = f"{input_bam.name}.pbi"
    num_reads = bam_utils.load_read_count(pbi) if os.path.exists(pbi) else None
    if not num_reads:
        num_reads = bam_utils.get_read_count_from_bam_index(input_bam)
    if num_reads:
        logger.info(f"About to correct umis in %d reads.", num_reads)

    # silence message about the .bai file not being found
    pysam.set_verbosity(0)

    logger.info(f"Writing UMI corrected reads to: {output_bam}")
    logger.info(f"Writing reads with uncorrectable UMIs to: {reject_bam}")

    # Split reads into groups by locus
    cached_file_names = [f"{input_bam.name}.locus2reads.pickle",
                         os.path.basename(f"{input_bam.name}.locus2reads.pickle")]

    if cache_read_loci and input_bam.seekable() and any([os.path.exists(p) for p in cached_file_names]):
        resolved_cache_name = [p for p in cached_file_names if os.path.exists(p)][0]
        logger.info(f"Loading locus2reads cache: {resolved_cache_name}")
        st = time.time()
        locus2reads = pickle.load(open(f"{resolved_cache_name}", 'rb'))
        et = time.time()
        logger.info(f"done. (elapsed: {et - st}s)")

    else:
        logger.info("Creating locus -> read map...")
        locus2reads = extract_read_groups(input_bam, umi_length, pre_extracted, max_ccs_umi_length_delta,
                                          max_clr_umi_length_delta, umi_tag, gene_tag, eq_class_tag, num_reads)
        if cache_read_loci:
            # Cache the locus2reads data:
            resolved_cache_name = os.path.basename(f"{input_bam.name}.locus2reads.pickle")
            logger.info(f"Writing locus2reads cache: {resolved_cache_name}")
            st = time.time()
            pickle.dump(locus2reads, open(f"{resolved_cache_name}", 'wb'))
            et = time.time()
            logger.info(f"done. (elapsed: {et - st}s)")

    logger.info("Number of loci: %d", len(locus2reads))

    # Reset our position to the start of the input file for our second traversal:
    input_bam.seek(0)

    # Correct reads at each locus
    logger.info("Correcting reads at each locus...")
    read2umi = {}
    for locus in tqdm(locus2reads, desc="Processing each locus", unit=" locus", colour="green", position=1,
                      file=sys.stderr, total=len(locus2reads), leave=False, disable=not sys.stdin.isatty()):
        process_reads_at_locus(locus2reads[locus], read2umi, umi_length,
                               max_ccs_edit_dist, max_clr_edit_dist,
                               max_ccs_length_diff, max_clr_length_diff,
                               max_ccs_gc_diff, max_clr_gc_diff)

    num_corrected = 0
    num_rejected = 0
    total_reads = 0

    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as input_bam_file:
        with pysam.AlignmentFile(output_bam, "wb", template=input_bam_file) as correct_umi_bam:
            with pysam.AlignmentFile(reject_bam, "wb", template=input_bam_file) as rejected_out_umi_bam:

                # Output BAM with corrected UMIs
                for read in tqdm(input_bam_file, desc="Writing out UMI-corrected reads", unit=" read", colour="green",
                                 position=2, file=sys.stderr, total=num_reads, leave=False,
                                 disable=not sys.stdin.isatty()):
                    if read.qname in read2umi:
                        read.set_tag(final_umi_tag, read2umi[read.qname])
                        read.set_tag(umi_corrected_tag, 1)
                    else:
                        read.set_tag(final_umi_tag, read.get_tag(umi_tag))
                        read.set_tag(umi_corrected_tag, 0)

                    if read_passes_filters(read, umi_length, min_back_seg_score, max_final_ccs_umi_length_delta,
                                           max_final_clr_umi_length_delta, final_umi_tag, back_alignment_score_tag):
                        correct_umi_bam.write(read)
                        num_corrected += 1
                    else:
                        rejected_out_umi_bam.write(read)
                        num_rejected += 1

                    total_reads += 1

    t_end = time.time()
    logger.info(f"Done. Elapsed time: {t_end - t_start:2.2f}s. "
                f"Overall processing rate: {total_reads / (t_end - t_start):2.2f} reads/s.")
    logger.info(f"Total Number of reads: {total_reads}")

    count_str, pct_str = get_field_count_and_percent_string(num_corrected, total_reads, FFORMAT)
    logger.info(f"Number of reads with corrected UMIs: {count_str} {pct_str}")

    count_str, pct_str = get_field_count_and_percent_string(num_rejected, total_reads, FFORMAT)
    logger.info(f"Number of reads with uncorrectable UMIs: {count_str} {pct_str}")


# ========================================================


def get_read_type(read):
    return ReadType.CCS if read.get_tag(READ_QUALITY_TAG) != -1 else ReadType.CLR


def get_read_locus(read, eq_class_tag):
    return read.get_tag(longbow.utils.constants.READ_BARCODE_CORRECTED_TAG), read.get_tag(eq_class_tag)


def get_read_seq(read, pre_extracted):
    if pre_extracted:
        return read.query_sequence.upper()
    else:
        seg = SegmentInfo.from_tag(read.get_tag(longbow.utils.constants.SEGMENTS_TAG))
        return read.query_sequence.upper()[seg.start:seg.end + 1]


def get_back_aln_score(read, back_alignment_score_tag):
    return int(read.get_tag(back_alignment_score_tag).split("/")[0])


def valid_umi(read, umi_length, ccs_max_umi_len_delta, clr_max_umi_len_delta, umi_tag):
    # checks the deviation of the UMI length
    if ReadType(get_read_type(read)) == ReadType.CCS:
        return abs(len(read.get_tag(umi_tag)) - umi_length) <= ccs_max_umi_len_delta
    else:
        return abs(len(read.get_tag(umi_tag)) - umi_length) <= clr_max_umi_len_delta


def valid_gene(read, gene_tag):
    # requires either a MAS-seq or Gencode gene tag
    return (MAS_GENE_PREFIX in read.get_tag(gene_tag)) or (GENCODE_GENE_PREFIX in read.get_tag(gene_tag))


def valid_tags(read, umi_tag, eq_class_tag):
    # checks for the presence of required tags
    return read.has_tag(READ_QUALITY_TAG) and read.has_tag(longbow.utils.constants.READ_BARCODE_CORRECTED_TAG) \
           and read.has_tag(umi_tag) and read.has_tag(eq_class_tag)


def read_passes_filters(read, umi_length, min_back_seg_score, max_final_ccs_umi_length_delta,
                        max_final_clr_umi_length_delta, final_umi_tag, back_alignment_score_tag):
    # filters the read based on the final UMI length and back alignment score

    max_umi_delta_filter = max_final_ccs_umi_length_delta \
        if ReadType(get_read_type(read)) == ReadType.CCS else max_final_clr_umi_length_delta
    return get_back_aln_score(read, back_alignment_score_tag) >= min_back_seg_score and abs(len(read.get_tag(final_umi_tag)) - umi_length) <= max_umi_delta_filter


def extract_read_groups(input_bam_fname, umi_length, pre_extracted, ccs_max_umi_len_delta, clr_max_umi_len_delta, umi_tag,
                        gene_tag, eq_class_tag, total_num_reads=None):
    locus2reads = defaultdict(list)
    n_filtered_umi = 0
    n_filtered_gene = 0
    n_valid_reads = 0
    n_total_reads = 0
    with pysam.AlignmentFile(input_bam_fname, "rb") as input_bam:
        for read in tqdm(input_bam, desc="Extracting Read Groups", unit=" read", total=total_num_reads, position=0):
            if not valid_tags(read, umi_tag, eq_class_tag):
                continue
            if not valid_gene(read, gene_tag):
                n_filtered_gene += 1
                continue
            if not valid_umi(read, umi_length, ccs_max_umi_len_delta, clr_max_umi_len_delta, umi_tag):
                n_filtered_umi += 1
                continue
            n_valid_reads += 1
            locus = get_read_locus(read, eq_class_tag)
            locus2reads[locus].append(ReadSnapshot(read, pre_extracted, umi_tag))
        logger.info("Number of reads: %d", n_total_reads)
        logger.info("Number of valid reads: %d", n_valid_reads)
        logger.info("Number of filtered by gene: %d", n_filtered_gene)
        logger.info("Number of filtered by UMI: %d", n_filtered_umi)
    return locus2reads


def get_conversion_type(source, target):
    return ReadType.CCS if ((source.type == target.type) and (source.type == ReadType.CCS)) else ReadType.CLR


def can_convert(source, target, max_ccs_edit_dist, max_clr_edit_dist, max_ccs_length_diff, max_clr_length_diff,
                max_ccs_gc_diff, max_clr_gc_diff):
    conversion_type = get_conversion_type(source, target)
    max_edit_dist = max_ccs_edit_dist if conversion_type == ReadType.CCS else max_clr_edit_dist

    edit_dist = levenshtein(source.umi, target.umi, max_edit_dist)
    if edit_dist > max_edit_dist:
        return False

    delta_len = abs(source.len - target.len)
    delta_gc = abs(source.gc - target.gc)
    max_len_diff = max_ccs_length_diff if conversion_type == ReadType.CCS else max_clr_length_diff
    max_gc_diff = max_ccs_gc_diff if conversion_type == ReadType.CCS else max_clr_gc_diff

    return (delta_len <= max_len_diff) and (delta_gc <= max_gc_diff)


def get_unique_targets(reads):
    return list(set(reads))


def build_graph(reads, max_ccs_edit_dist, max_clr_edit_dist, max_ccs_length_diff, max_clr_length_diff, max_ccs_gc_diff,
                max_clr_gc_diff):
    targets = get_unique_targets(reads)
    graph = defaultdict(list)
    target2umi_counts = defaultdict(Counter)
    target2umi_seq = {target_id: target.umi for target_id, target in enumerate(targets)}

    # For very large groups we need to tell the user that somehting is actually happening:
    if len(reads) * len(targets) < 1e6:
        for read_id, read in enumerate(reads):
            for target_id, target in enumerate(targets):
                if can_convert(read, target, max_ccs_edit_dist, max_clr_edit_dist, max_ccs_length_diff,
                               max_clr_length_diff, max_ccs_gc_diff, max_clr_gc_diff):
                    graph[target_id].append(read_id)
                    target2umi_counts[target_id][read.umi] += 1
    else:
        for read_id, read in tqdm(enumerate(reads), desc="Building graph", position=2):
            for target_id, target in enumerate(targets):
                if can_convert(read, target, max_ccs_edit_dist, max_clr_edit_dist, max_ccs_length_diff,
                               max_clr_length_diff, max_ccs_gc_diff, max_clr_gc_diff):
                    graph[target_id].append(read_id)
                    target2umi_counts[target_id][read.umi] += 1

    return graph, target2umi_counts, target2umi_seq


def min_vertex_cover(read_ids, graph, target2umi_counts, target2umi_seq):
    umi_groups = []
    while len(read_ids) > 0:
        # find the largest group, tie breaking: target matching the max support UMI in group
        max_target_id = max(graph.keys(),
                            key=lambda t: (len(graph[t]), max(target2umi_counts[t],
                                                              key=target2umi_counts[t].get) == target2umi_seq[t]))
        max_size = len(graph[max_target_id])

        if max_size == 1:
            break

        umi_groups.append(graph[max_target_id])

        # remove reads in the largest group from other groups
        for selected_read_id in graph[max_target_id]:
            for target_id in graph:
                if target_id == max_target_id:
                    continue
                for read_id in graph[target_id]:
                    if selected_read_id == read_id:
                        graph[target_id].remove(read_id)
            read_ids.remove(selected_read_id)

        del graph[max_target_id]

    return umi_groups


def process_reads_at_locus(reads, read2umi, umi_length,
                           max_ccs_edit_dist, max_clr_edit_dist,
                           max_ccs_length_diff, max_clr_length_diff,
                           max_ccs_gc_diff, max_clr_gc_diff):
    if len(reads) < 2:
        return

    # if all the reads have the same umi, skip correction
    if all(read.umi == reads[0].umi for read in reads):
        return

    graph, target2umi_counts, target2umi_seq = build_graph(reads, max_ccs_edit_dist, max_clr_edit_dist,
                                                           max_ccs_length_diff, max_clr_length_diff, max_ccs_gc_diff,
                                                           max_clr_gc_diff)

    read_ids = list(range(len(reads)))
    umi_groups = min_vertex_cover(read_ids, graph, target2umi_counts, target2umi_seq)
    for group in umi_groups:
        umis = [reads[read_id].umi for read_id in group]
        # pick a umi with maximal support and closest len to UMI_LEN
        umi_max = max(umis, key=lambda t: (umis.count(t), -abs(umi_length - len(t))))
        for read_id in group:
            read2umi[reads[read_id].name] = umi_max
