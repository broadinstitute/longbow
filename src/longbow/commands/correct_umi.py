import logging
import os
import re
import time
import sys
import itertools

import click
import click_log
from tqdm import tqdm
import enum
import operator

from polyleven import levenshtein

import pysam

from construct import *

from collections import defaultdict
from collections import Counter

import longbow.utils.constants
from ..utils import bam_utils

PROG_NAME = "correct_umi"

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(PROG_NAME)
click_log.basic_config(logger)

# Settings for the correction.
# These have been tuned by examining data.
MAX_EDIT_DISTANCE = {"CCS": 2, "CLR": 3}
MAX_LENGTH_DIFF = {"CCS": 50, "CLR": 150}
MAX_GC_CONTENT_DIFF = {"CCS": 0.05, "CLR": 0.15}
CONFIG_FILTER_OP = "AND"
MAX_UMI_DELTA = {"CCS": 3, "CLR": 4}
MAX_UMI_DELTA_FILTER = {"CCS": 3, "CLR": 3}
MIN_BACK_ALIGNMENT_SCORE = 10

MAS_GENE_PREFIX = "MAS"
GENCODE_GENE_PREFIX = "ENSG"


class ReadType(enum.Enum):
    CCS = enum.auto()
    CLR = enum.auto()


UMI_TAG = "JX"  # "ZU"
FINAL_UMI_TAG = "BX"
UMI_CORR_TAG = "UX"
EQ_CLASS_TAG = "eq"
GENE_TAG = "XG"
CODING_REGION = "cDNA"
READ_QUALITY_TAG = "rq"
BACK_ALIGNMENT_SCORE_TAG = "JB"


class ReadSnapshot:
    def __init__(self, read, pre_extracted):
        self.umi = read.get_tag(UMI_TAG)
        self.type = get_read_type(read)
        self.start = read.reference_start
        self.end = read.reference_end
        sequence = get_read_seq(read, pre_extracted)
        self.len = len(sequence)
        self.gc = float(sequence.count('C') + sequence.count('G'))/len(sequence)
        self.name = read.qname

    def __eq__(self, other):
        return self.umi == other.umi and self.type == other.type and self.start == other.start and \
               self.end == other.end and self.len == other.len and \
               self.gc == other.gc and self.name == other .name


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
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
def main(umi_length, output_bam, reject_bam, force, pre_extracted, input_bam):
    """Correct UMIs with Set Cover algorithm."""

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

    # split reads into groups by locus
    logger.info("Creating locus -> read map...")
    locus2reads = create_read_loci(input_bam, umi_length, pre_extracted)
    logger.info("Number of loci: %d", len(locus2reads))

    # Reset our position to the start of the input file for our second traversal:
    input_bam.seek(0)

    # Correct reads at each locus
    logger.info("Correcting reads at each locus...")
    read2umi = {}
    for locus in tqdm(locus2reads, desc="Processing each locus", unit=" locus", colour="green",
                      file=sys.stderr, total=len(locus2reads), leave=False, disable=not sys.stdin.isatty()):
        process_reads_at_locus(locus2reads[locus], read2umi, umi_length)

    num_corrected = 0
    num_rejected = 0
    total_reads = 0

    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as input_bam_file:
        with pysam.AlignmentFile(output_bam, "wb", template=input_bam_file) as correct_umi_bam:
            with pysam.AlignmentFile(reject_bam, "wb", template=input_bam_file) as rejected_out_umi_bam:

                # Output BAM with corrected UMIs
                for read in tqdm(input_bam_file, desc="Writing out UMI-corrected reads", unit=" read", colour="green",
                                 file=sys.stderr, total=num_reads, leave=False, disable=not sys.stdin.isatty()):
                    if read.qname in read2umi:
                        read.set_tag(FINAL_UMI_TAG, read2umi[read.qname])
                        read.set_tag(UMI_CORR_TAG, 1)
                    else:
                        read.set_tag(FINAL_UMI_TAG, read.get_tag(UMI_TAG))
                        read.set_tag(UMI_CORR_TAG, 0)

                    if read_passes_filters(read, umi_length):
                        correct_umi_bam.write(read)
                        num_corrected += 1
                    else:
                        rejected_out_umi_bam.write(read)
                        num_rejected += 1

                    total_reads += 1

    t_end = time.time()
    logger.info(f"Done. Elapsed time: {t_end - t_start:2.2f}s. "
                f"Overall processing rate: {total_reads / (t_end - t_start):2.2f} reads/s.")
    logger.info(f"Number of reads with corrected UMIs: {num_corrected}/{total_reads} "
                f"({100*(num_corrected/total_reads):2.2f}%)")
    logger.info(f"Number of reads with uncorrectable UMIs: {num_rejected}/{total_reads} "
                f"({100 * (num_rejected / total_reads):2.2f}%)")
    logger.info(f"Total Number of reads: {total_reads}")


# ========================================================


def get_read_type(read):
    return ReadType.CCS if read.get_tag(READ_QUALITY_TAG) != -1 else ReadType.CLR


def get_read_locus(read):
    return read.get_tag(longbow.utils.constants.READ_BARCODE_CORRECTED_TAG), read.get_tag(EQ_CLASS_TAG)


def get_read_seq(read, pre_extracted):
    if pre_extracted:
        return read.query_sequence.upper()
    else:
        start, end = read.get_tag(longbow.utils.constants.SEGMENTS_TAG).split(f"{CODING_REGION}:", 1)[1].split(",")[
            0].split("-")
        return read.query_sequence.upper()[int(start):int(end) + 1]


def get_back_aln_score(read):
    return int(read.get_tag(BACK_ALIGNMENT_SCORE_TAG).split("/")[0])


def valid_umi(read, umi_length):
    # checks the deviation of the UMI length
    return abs(len(read.get_tag(UMI_TAG)) - umi_length) <= MAX_UMI_DELTA[ReadType(get_read_type(read)).name]


def valid_gene(read):
    # requires either a MAS-seq or Gencode gene tag
    return (MAS_GENE_PREFIX in read.get_tag(GENE_TAG)) or (GENCODE_GENE_PREFIX in read.get_tag(GENE_TAG))


def valid_tags(read):
    # checks for the presence of required tags
    return read.has_tag(READ_QUALITY_TAG) and read.has_tag(longbow.utils.constants.READ_BARCODE_CORRECTED_TAG) \
           and read.has_tag(UMI_TAG) and read.has_tag(EQ_CLASS_TAG)


def read_passes_filters(read, umi_length):
    # filters the read based on the final UMI length and back alignment score
    return get_back_aln_score(read) >= MIN_BACK_ALIGNMENT_SCORE and \
           abs(len(read.get_tag(FINAL_UMI_TAG)) - umi_length) <= MAX_UMI_DELTA_FILTER[ReadType(get_read_type(read)).name]


def create_read_loci(input_bam_fname, umi_length, pre_extracted):
    locus2reads = defaultdict(list)
    n_filtered_umi = 0
    n_filtered_gene = 0
    n_valid_reads = 0
    with pysam.AlignmentFile(input_bam_fname, "rb") as input_bam:
        for read in tqdm(input_bam, desc="Extracting Read Groups", unit=" read"):
            if not valid_tags(read):
                continue
            if not valid_gene(read):
                n_filtered_gene += 1
                continue
            if not valid_umi(read, umi_length):
                n_filtered_umi += 1
                continue
            n_valid_reads += 1
            locus = get_read_locus(read)
            locus2reads[locus].append(ReadSnapshot(read, pre_extracted))
        logger.info("Number of valid reads: %d", n_valid_reads)
        logger.info("Number of filtered by gene: %d", n_filtered_gene)
        logger.info("Number of filtered by UMI: %d", n_filtered_umi)
    return locus2reads


def get_conversion_type(source, target):
    return "CCS" if ((source.type == target.type) and (source.type == ReadType.CCS)) else "CLR"


def can_convert(source, target):
    conversion_type = get_conversion_type(source, target)
    edit_dist = levenshtein(source.umi, target.umi, MAX_EDIT_DISTANCE[conversion_type])
    if edit_dist > MAX_EDIT_DISTANCE[conversion_type]:
        return False
    delta_len = abs(source.len - target.len)
    delta_gc = abs(source.gc - target.gc)
    op = operator.and_
    if op == "OR":
        op = operator.or_
    return op(delta_len <= MAX_LENGTH_DIFF[conversion_type],
              delta_gc <= MAX_GC_CONTENT_DIFF[conversion_type])


def get_unique_targets(reads):
    return list(set(reads))


def build_graph(reads):
    targets = get_unique_targets(reads)
    graph = defaultdict(list)
    target2umi_counts = defaultdict(Counter)
    target2umi_seq = {target_id: target.umi for target_id, target in enumerate(targets)}
    for read_id, read in enumerate(reads):
        for target_id, target in enumerate(targets):
            if can_convert(read, target):
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


def process_reads_at_locus(reads, read2umi, umi_length):
    if len(reads) < 2:
        return

    # if all the reads have the same umi, skip correction
    if all(read.umi == reads[0].umi for read in reads):
        return

    graph, target2umi_counts, target2umi_seq = build_graph(reads)
    read_ids = list(range(len(reads)))
    umi_groups = min_vertex_cover(read_ids, graph, target2umi_counts, target2umi_seq)
    for group in umi_groups:
        umis = [reads[read_id].umi for read_id in group]
        # pick a umi with maximal support and closest len to UMI_LEN
        umi_max = max(umis, key=lambda t: (umis.count(t), -abs(umi_length - len(t))))
        for read_id in group:
            read2umi[reads[read_id].name] = umi_max
