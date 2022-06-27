import logging
import queue
import time
import sys
import os
import itertools
import tempfile
import enum

import click
import click_log

import pysam
import multiprocessing as mp

from collections import defaultdict

from tqdm import tqdm

from construct import *

import longbow.utils.constants
from ..utils import bam_utils, barcode_utils
from ..utils.model import LibraryModel


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("correct")
click_log.basic_config(logger)


class BarcodeResolutionFailure(enum.Enum):
    NO_RAW_BARCODE = enum.auto()
    AMBIGUOUS = enum.auto()
    NO_MATCH_IN_LEV_DIST = enum.auto()


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-p",
    "--pbi",
    required=False,
    type=click.Path(),
    help="BAM .pbi index file",
)
@click.option(
    "-t",
    "--threads",
    type=int,
    default=mp.cpu_count() - 1,
    show_default=True,
    help="number of threads to use (0 for all)",
)
@click.option(
    "-o",
    "--output-bam",
    default="-",
    type=click.Path(exists=False),
    help="annotated bam output  [default: stdout]",
)
@click.option(
    "-m",
    "--model",
    help="The model to use for annotation.  If not specified, it will be autodetected from "
         "the BAM header.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name "
         "and Longbow will attempt to read in the file and create a LibraryModel from it.  "
         "Longbow will assume the contents are the configuration of a LibraryModel as per "
         "LibraryModel.to_json()."
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.option(
    '-r',
    '--restrict-to-allowlist',
    is_flag=True,
    default=True,
    show_default=True,
    help="Restrict barcode correction possibilities to only those on the allowlist."
)
@click.option(
    "-b",
    "--barcode-tag",
    type=str,
    default=longbow.utils.constants.READ_BARCODE_TAG,
    show_default=True,
    help="The tag from which to read the uncorrected barcode."
)
@click.option(
    "-c",
    "--corrected-tag",
    type=str,
    default=longbow.utils.constants.READ_BARCODE_CORRECTED_TAG,
    show_default=True,
    help="The tag in which to store the corrected barcode."
)
@click.option(
    "-a",
    "--allow-list",
    required=True,
    type=click.Path(exists=True),
    help="List of allowed barcodes for specified tag (.txt, .txt.gz).",
)
@click.option(
    "--barcode-freqs",
    type=click.Path(exists=True),
    required=False,
    help="TSV file containing barcodes and the frequencies associated with them in the data (BARCODE\tFREQ).  "
         "If not provided, barcode freqs will be uniformly seeded by the barcode whitelist.  "
         "NOTE: If barcodes in this freqs file are not in the allow list and the `-r` flag is not given, it is possible"
         "to end up with reads that have barcodes which where corrected to values that are not on the allow list.",
)
@click.option(
    "--max-hifi-dist",
    type=int,
    default=2,
    show_default=True,
    help="Maximum levenshtein distance to allow for hifi/CCS reads during correction.",
)
@click.option(
    "--max-clr-dist",
    type=int,
    default=3,
    show_default=True,
    help="Maximum levenshtein distance to allow for CLR (ccs uncorrected) reads during correction.",
)
@click.option(
    "--ccs-corrected-rq-threshold",
    type=float,
    default=0,
    show_default=True,
    help="Value of the `rq` tag above which reads are considered to be CCS corrected (Hifi).",
)
@click.option(
    "--barcode-uncorrectable-bam",
    default="/dev/null",
    show_default=True,
    type=click.Path(exists=False),
    help="File to which to write all reads with barcodes that could not be corrected.",
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, threads, output_bam, model, force, restrict_to_allowlist, barcode_tag, corrected_tag, allow_list,
         barcode_freqs, max_hifi_dist, max_clr_dist, ccs_corrected_rq_threshold, barcode_uncorrectable_bam, input_bam):
    """Correct tag to values provided in barcode allowlist."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    # Decide if we're disabling progress bars:
    disable_pbar = not sys.stdin.isatty()

    logger.info(f"Writing reads with corrected barcodes to: {output_bam}")
    logger.info(f"Writing reads with barcodes that could not be corrected to: {barcode_uncorrectable_bam}")

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    logger.info(f"Using CCS Levenshtein distance threshold: {max_hifi_dist}")
    logger.info(f"Using CLR Levenshtein distance threshold: {max_clr_dist}")

    # Load barcode allow list:
    barcode_allow_list = barcode_utils.load_barcode_allowlist(allow_list, disable_pbar=disable_pbar)

    # Load number of reads, if pbi exists:
    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    num_reads = bam_utils.load_read_count(pbi) if os.path.exists(pbi) else None

    # Create queues for data:
    # NOTE: We're using processes to overcome the Global Interpreter Lock.
    queue_size = threads * 2 if threads < 10 else 20
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(maxsize=queue_size)
    # process_input_data_queue = manager.Queue()
    results = manager.Queue()

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file:

        # Load our model:
        if model is None:
            lb_model = LibraryModel.from_json_obj(bam_utils.get_model_from_bam_header(bam_file.header))
        elif model is not None and LibraryModel.has_prebuilt_model(model):
            lb_model = LibraryModel.build_pre_configured_model(model)
        else:
            lb_model = LibraryModel.from_json_file(model)
        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        # Get our barcode length:
        barcode_length = _get_barcode_tag_length_from_model(lb_model, barcode_tag)

        # Generate our symspellpy index:
        if not barcode_freqs:
            logger.info("No barcode freq file provided.")
            barcode_freqs = _generate_barcode_freqs_from_allow_list(barcode_allow_list)
            logger.info(f"Generated uniform frequencies from allow list in file: {barcode_freqs}")
        else:
            logger.info(f"Using barcode frequencies from: {barcode_freqs}")
            if restrict_to_allowlist:
                logger.info("Filtering barcode frequencies file to only those on the allow list...")
                st = time.time()
                barcode_freqs = _filter_barcode_freqs_by_allow_list(barcode_freqs, barcode_allow_list)
                logger.info(f"Barcode freqs filtered in {time.time() - st:2.4f}s")
            else:
                logger.warning("Allowing ANY barcode in freqs file to seed corrections.  "
                               "THIS MAY PRODUCE NON-ALLOWLIST BARCODES IN RESULTING CORRECTED DATA!")

        # Start worker sub-processes:
        worker_process_pool = []
        for _ in range(threads):
            p = mp.Process(
                target=_correct_barcode_fn,
                args=(process_input_data_queue, results, bam_file.header.to_dict(), barcode_tag, corrected_tag,
                      barcode_length, barcode_allow_list, barcode_freqs, max_hifi_dist, max_clr_dist,
                      ccs_corrected_rq_threshold)
            )
            p.start()
            worker_process_pool.append(p)

        # Get header for output file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[lb_model])

        # Start output worker:
        res = manager.dict({"num_ccs_reads": 0, "num_ccs_reads_corrected": 0, "num_ccs_reads_raw_was_correct": 0,
                            "num_clr_reads": 0, "num_clr_reads_corrected": 0, "num_clr_reads_raw_was_correct": 0,
                            "num_ccs_with_barcodes": 0, "num_clr_with_barcodes": 0,
                            "num_ccs_could_not_correct_ambiguous": 0, "num_ccs_could_not_correct_no_match": 0,
                            "num_clr_could_not_correct_ambiguous": 0, "num_clr_could_not_correct_no_match": 0, })
        output_worker = mp.Process(
            target=_write_thread_fn,
            args=(results, out_header, output_bam, barcode_uncorrectable_bam, barcode_tag, corrected_tag,
                  ccs_corrected_rq_threshold, num_reads, disable_pbar, res),
        )
        output_worker.start()

        # Add in a sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        iter_data = itertools.chain(bam_file, (None,) * threads)
        for r in iter_data:
            if r is not None:
                process_input_data_queue.put(r.to_string())
            else:
                process_input_data_queue.put(r)

        # Wait for our input jobs to finish:
        for p in worker_process_pool:
            p.join()

        # Now that our input processes are done, we can add our exit sentinel onto the output queue and
        # wait for that process to end:
        results.put(None)
        output_worker.join()

    # Print out some stats:

    stat_prefix = "STATS: "

    total_reads = res["num_ccs_reads"] + res["num_clr_reads"]
    logger.info(f"{stat_prefix}Total reads seen: {total_reads}")

    total_with_barcodes = res["num_ccs_with_barcodes"] + res["num_clr_with_barcodes"]
    count_str, pct_str = _get_field_count_and_percent_string(total_with_barcodes, total_reads)
    logger.info(f"{stat_prefix}Reads with {barcode_tag} barcodes: {count_str} {pct_str}")

    total_without_barcodes = total_reads - total_with_barcodes
    count_str, pct_str = _get_field_count_and_percent_string(total_without_barcodes, total_reads)
    logger.info(f"{stat_prefix}Reads without {barcode_tag} barcodes: {count_str} {pct_str}")

    total_corrected = res["num_ccs_reads_corrected"] + res["num_clr_reads_corrected"]
    count_str, pct_str = _get_field_count_and_percent_string(total_corrected, total_reads)
    logger.info(f"{stat_prefix}Reads able to be corrected into {corrected_tag} tag: {count_str} {pct_str}")

    total_uncorrected = total_reads - (res["num_ccs_reads_corrected"] + res["num_clr_reads_corrected"])
    count_str, pct_str = _get_field_count_and_percent_string(total_uncorrected, total_reads)
    logger.info(f"{stat_prefix}Reads unable to be corrected into {corrected_tag} tag: {count_str} {pct_str}")

    total_uncorrected_ambiguous = res["num_ccs_could_not_correct_ambiguous"] + res["num_clr_could_not_correct_ambiguous"]
    count_str, pct_str = _get_field_count_and_percent_string(total_uncorrected_ambiguous, total_reads)
    count_str2, pct_str2 = _get_field_count_and_percent_string(total_uncorrected_ambiguous, total_uncorrected)
    logger.info(f"{stat_prefix}Reads unable to be corrected into {corrected_tag} tag - Ambiguous (of total reads): {count_str} {pct_str}")
    logger.info(f"{stat_prefix}Reads unable to be corrected into {corrected_tag} tag - Ambiguous (of uncorrected reads): {count_str2} {pct_str2}")

    total_uncorrected_no_match = res["num_ccs_could_not_correct_no_match"] + res["num_clr_could_not_correct_no_match"]
    count_str, pct_str = _get_field_count_and_percent_string(total_uncorrected_no_match, total_reads)
    count_str2, pct_str2 = _get_field_count_and_percent_string(total_uncorrected_no_match, total_uncorrected)
    logger.info(f"{stat_prefix}Reads unable to be corrected into {corrected_tag} tag - No Match in Lev Dist (of total reads): {count_str} {pct_str}")
    logger.info(f"{stat_prefix}Reads unable to be corrected into {corrected_tag} tag - No Match in Lev Dist (of uncorrected reads): {count_str2} {pct_str2}")

    total_were_already_correct = res["num_ccs_reads_raw_was_correct"] + res["num_clr_reads_raw_was_correct"]
    count_str, pct_str = _get_field_count_and_percent_string(total_were_already_correct, total_reads)
    logger.info(f"{stat_prefix}Reads with already correct {barcode_tag} barcodes: {count_str} {pct_str}")

    for read_type in ['ccs', 'clr']:

        logger.info("=" * 80)

        count_str, pct_str = _get_field_count_and_percent_string(res[f"num_{read_type}_reads"], total_reads)
        logger.info(f"{stat_prefix}{read_type.upper()} reads seen: {count_str} {pct_str}")

        count_str, pct_str = _get_field_count_and_percent_string(res[f"num_{read_type}_with_barcodes"], total_reads)
        count_str2, pct_str2 = _get_field_count_and_percent_string(res[f"num_{read_type}_with_barcodes"], res[f"num_{read_type}_reads"])
        logger.info(f"{stat_prefix}{read_type.upper()} reads with {barcode_tag} barcodes (of total reads): {count_str} {pct_str}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads with {barcode_tag} barcodes (of {read_type} reads): {count_str2} {pct_str2}")

        count_str, pct_str = _get_field_count_and_percent_string(res[f"num_{read_type}_reads"] - res[f"num_{read_type}_with_barcodes"], total_reads)
        count_str2, pct_str2 = _get_field_count_and_percent_string(res[f"num_{read_type}_reads"] - res[f"num_{read_type}_with_barcodes"], res[f"num_{read_type}_reads"])
        logger.info(f"{stat_prefix}{read_type.upper()} reads without {barcode_tag} barcodes (of total reads): {count_str} {pct_str}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads without {barcode_tag} barcodes (of {read_type} reads): {count_str2} {pct_str2}")

        count_str, pct_str = _get_field_count_and_percent_string(res[f"num_{read_type}_reads_corrected"], total_reads)
        count_str2, pct_str2 = _get_field_count_and_percent_string(res[f"num_{read_type}_reads_corrected"], res[f"num_{read_type}_reads"])
        logger.info(f"{stat_prefix}{read_type.upper()} reads able to be corrected into {corrected_tag} tag (of total reads): {count_str} {pct_str}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads able to be corrected into {corrected_tag} tag (of {read_type} reads): {count_str2} {pct_str2}")

        count_str, pct_str = _get_field_count_and_percent_string(res[f"num_{read_type}_reads"] - res[f"num_{read_type}_reads_corrected"], total_reads)
        count_str2, pct_str2 = _get_field_count_and_percent_string(res[f"num_{read_type}_reads"] - res[f"num_{read_type}_reads_corrected"], res[f"num_{read_type}_reads"])
        logger.info(f"{stat_prefix}{read_type.upper()} reads unable to be corrected into {corrected_tag} tag (of total reads): {count_str} {pct_str}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads unable to be corrected into {corrected_tag} tag (of {read_type} reads): {count_str2} {pct_str2}")

        count_str, pct_str = _get_field_count_and_percent_string(res[f"num_{read_type}_could_not_correct_ambiguous"], total_reads)
        count_str2, pct_str2 = _get_field_count_and_percent_string(res[f"num_{read_type}_could_not_correct_ambiguous"], res[f"num_{read_type}_reads"])
        count_str3, pct_str3 = _get_field_count_and_percent_string(res[f"num_{read_type}_could_not_correct_ambiguous"], res[f"num_{read_type}_reads"] - res[f"num_{read_type}_reads_corrected"])
        logger.info(f"{stat_prefix}{read_type.upper()} reads unable to be corrected into {corrected_tag} tag - Ambiguous (of total reads): {count_str} {pct_str}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads unable to be corrected into {corrected_tag} tag - Ambiguous (of {read_type} reads): {count_str2} {pct_str2}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads unable to be corrected into {corrected_tag} tag - Ambiguous (of {read_type} uncorrected reads): {count_str3} {pct_str3}")

        count_str, pct_str = _get_field_count_and_percent_string(res[f"num_{read_type}_could_not_correct_no_match"], total_reads)
        count_str2, pct_str2 = _get_field_count_and_percent_string(res[f"num_{read_type}_could_not_correct_no_match"], res[f"num_{read_type}_reads"])
        count_str3, pct_str3 = _get_field_count_and_percent_string(res[f"num_{read_type}_could_not_correct_no_match"], res[f"num_{read_type}_reads"] - res[f"num_{read_type}_reads_corrected"])
        logger.info(f"{stat_prefix}{read_type.upper()} reads unable to be corrected into {corrected_tag} tag - No Match in Lev Dist (of total reads): {count_str} {pct_str}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads unable to be corrected into {corrected_tag} tag - No Match in Lev Dist (of {read_type} reads): {count_str2} {pct_str2}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads unable to be corrected into {corrected_tag} tag - No Match in Lev Dist (of {read_type} uncorrected reads): {count_str3} {pct_str3}")

        count_str, pct_str = _get_field_count_and_percent_string(res[f"num_{read_type}_reads_raw_was_correct"], total_reads)
        count_str2, pct_str2 = _get_field_count_and_percent_string(res[f"num_{read_type}_reads_raw_was_correct"], res[f"num_{read_type}_reads"])
        logger.info(f"{stat_prefix}{read_type.upper()} reads with already correct {corrected_tag} tag (of total reads): {count_str} {pct_str}")
        logger.info(f"{stat_prefix}{read_type.upper()} reads with already correct {corrected_tag} tag (of {read_type} reads): {count_str2} {pct_str2}")

    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {total_reads/(et - t_start):2.2f} reads/s.")


def _get_field_count_and_percent_string(count, total):
    count_str = f"{count}/{total}"
    pct_str = f"({0:2.4f}%)" if total == 0 else f"({100.0*count/total:2.4f}%)"

    return count_str, pct_str


def _write_thread_fn(data_queue, out_bam_header, out_bam_file_name, barcode_uncorrectable_bam,
                     barcode_tag, corrected_tag, ccs_corrected_rq_threshold, num_reads, disable_pbar, res):
    """Thread / process fn to write out all our data."""

    with pysam.AlignmentFile(out_bam_file_name, "wb", header=out_bam_header) as out_bam_file, \
        pysam.AlignmentFile(barcode_uncorrectable_bam, "wb", header=out_bam_header) as barcode_uncorrectable_bam, \
        tqdm(
            desc="Progress",
            unit=" read",
            colour="green",
            file=sys.stderr,
            disable=disable_pbar,
            total=num_reads,
            leave=False
        ) as pbar:

        while True:
            # Wait for some output data:
            raw_data = data_queue.get()

            # Check for exit sentinel:
            if raw_data is None:
                break
            # Should really never be None, but just in case:
            elif raw_data is None:
                continue

            # Unpack data:
            read, num_segments, num_corrected_segments, has_barcode = raw_data
            read = pysam.AlignedSegment.fromstring(read, out_bam_header)

            # Save ourselves some trouble later:
            ccs_type_stat_string = "ccs" if read.get_tag("rq") > ccs_corrected_rq_threshold else "clr"

            # Write our our read:
            if read.get_tag(longbow.utils.constants.COULD_CORRECT_BARCODE_TAG):
                out_bam_file.write(read)
            else:
                # "num_ccs_could_not_correct_ambiguous": 0, "num_ccs_could_not_correct_no_match"
                # Determine what the problem in correction was for our stats:
                correction_failure = read.get_tag(corrected_tag)

                # NOTE: using the enum.value field here because of the TYPE of the barcode tag.
                if correction_failure == BarcodeResolutionFailure.AMBIGUOUS.value:
                    res[f"num_{ccs_type_stat_string}_could_not_correct_ambiguous"] += 1
                elif correction_failure == BarcodeResolutionFailure.NO_MATCH_IN_LEV_DIST.value:
                    res[f"num_{ccs_type_stat_string}_could_not_correct_no_match"] += 1
                elif correction_failure == BarcodeResolutionFailure.NO_RAW_BARCODE.value:
                    # This case is already accounted for.
                    pass
                else:
                    raise RuntimeError(f"Unknown BarcodeResolutionFailure type ({correction_failure})!  "
                                       f"This should never happen!")

                # Set our corrected barcode tag to the unlabeled placeholder:
                read.set_tag(corrected_tag, longbow.utils.constants.UNLABELED_BARCODE)
                barcode_uncorrectable_bam.write(read)

            # Increment our counters:
            res[f"num_{ccs_type_stat_string}_reads"] += 1
            if read.has_tag(barcode_tag):
                res[f"num_{ccs_type_stat_string}_with_barcodes"] += 1
            if read.get_tag(longbow.utils.constants.COULD_CORRECT_BARCODE_TAG) == 1:
                res[f"num_{ccs_type_stat_string}_reads_corrected"] += 1
                if read.get_tag(longbow.utils.constants.BARCODE_CORRECTION_PERFORMED) == 0:
                    res[f"num_{ccs_type_stat_string}_reads_raw_was_correct"] += 1

            pbar.update(1)


def _correct_barcode_fn(in_queue, out_queue, bam_header_dict, barcode_tag, corrected_tag, barcode_length,
                        bc_allow_list, barcode_freqs, max_hifi_dist, max_clr_dist, ccs_corrected_rq_threshold):
    """Function to run in each subprocess.
    Replace barcode with corrected value."""

    bam_header = pysam.AlignmentHeader.from_dict(bam_header_dict)

    # Create the sym spell index here (in the sub-process) to reduce memory overhead / copying.
    logger.info(f"Generating barcode index...")
    st = time.time()
    sym_spell_index = barcode_utils.generate_symspell_index(barcode_freqs, max(max_hifi_dist, max_clr_dist),
                                                            barcode_length)
    logger.info(f"Barcode index generated in {time.time() - st:2.4f}s")

    while True:
        # Wait until we get some data.
        # Note: Because we have a sentinel value None inserted at the end of the input data for each
        #       subprocess, we don't have to add a timeout - we're guaranteed each process will always have
        #       at least one element.
        # try:
        #     raw_data = in_queue.get_nowait()
        # except queue.Empty:
        #     time.sleep(0.01)
        #     continue
        raw_data = in_queue.get()

        # Check for exit sentinel:
        if raw_data is None:
            return

        # Unpack our data here:
        read = pysam.AlignedSegment.fromstring(raw_data, bam_header)

        num_segments = 0
        num_corrected_segments = 0
        has_barcode = False

        if read.has_tag(barcode_tag):

            has_barcode = True
            old_bc = read.get_tag(barcode_tag)
            dist_threshold = max_hifi_dist if read.get_tag('rq') > ccs_corrected_rq_threshold else max_clr_dist

            offset = 0

            # Check to see if we have a barcode that is longer than expected.
            # If we do, then we should try perform matches at every valid position and take the best result:
            if len(old_bc) > barcode_length:
                new_bc, result_status, offset = _perform_barcode_multi_match(barcode_length, bc_allow_list,
                                                                             dist_threshold, old_bc, sym_spell_index)
            else:
                new_bc, _, result_status = \
                    barcode_utils.find_match_symspell(old_bc, bc_allow_list, sym_spell_index, dist_threshold)

            num_segments += 1
            if new_bc is not None:
                read.set_tag(corrected_tag, new_bc)
                read.set_tag(longbow.utils.constants.COULD_CORRECT_BARCODE_TAG, True)

                print(f"Old BC: {old_bc}")
                print(f"New BC: {new_bc}")

                # Need to do some basic math for padded barcode strings:
                if len(old_bc) > barcode_length:
                    end_pad_bases = int((len(old_bc) - barcode_length)/2)
                    correction_performed = new_bc != old_bc[end_pad_bases:end_pad_bases+barcode_length]
                else:
                    correction_performed = new_bc != old_bc
                read.set_tag(longbow.utils.constants.BARCODE_CORRECTION_PERFORMED, correction_performed)

                read.set_tag(longbow.utils.constants.READ_ADJUSTED_BARCODE_START, offset)

                num_corrected_segments += 1
            else:
                read.set_tag(longbow.utils.constants.COULD_CORRECT_BARCODE_TAG, False)
                read.set_tag(longbow.utils.constants.BARCODE_CORRECTION_PERFORMED, False)

                if result_status == barcode_utils.SymSpellMatchResultType.NO_MATCH_IN_LEV_DIST:
                    # NOTE: using the enum.value field here because of the TYPE of the barcode tag.
                    read.set_tag(corrected_tag, BarcodeResolutionFailure.NO_MATCH_IN_LEV_DIST.value)
                elif result_status == barcode_utils.SymSpellMatchResultType.AMBIGUOUS:
                    # NOTE: using the enum.value field here because of the TYPE of the barcode tag.
                    read.set_tag(corrected_tag, BarcodeResolutionFailure.AMBIGUOUS.value)
                else:
                    raise RuntimeError(f"Unknown SymSpellMatchResultType ({result_status})!  This should never happen!")
        else:
            read.set_tag(longbow.utils.constants.COULD_CORRECT_BARCODE_TAG, False)
            read.set_tag(longbow.utils.constants.BARCODE_CORRECTION_PERFORMED, False)
            # NOTE: using the enum.value field here because of the TYPE of the barcode tag.
            read.set_tag(corrected_tag, BarcodeResolutionFailure.NO_RAW_BARCODE.value)

        # Process and place our data on the output queue:
        out_queue.put(tuple([read.to_string(), num_segments, num_corrected_segments, has_barcode]))


def _perform_barcode_multi_match(barcode_length, bc_allow_list, dist_threshold, old_bc, sym_spell_index):
    # OK, we have to do multiple matches to try to find the barcode here:
    success_dict = dict()
    failure_counts = defaultdict(int)
    for i in range(0, len(old_bc) - barcode_length + 1):
        bc_sub_string = old_bc[i:i + barcode_length]

        new_bc, edit_dist, result_status = \
            barcode_utils.find_match_symspell(bc_sub_string, bc_allow_list, sym_spell_index, dist_threshold)
        if new_bc is not None:
            try:
                success_dict[edit_dist].append((new_bc, result_status, i))
            except KeyError:
                success_dict[edit_dist] = [(new_bc, result_status, i)]
        else:
            failure_counts[result_status] += 1

    # Setup for later:
    new_bc = None
    result_status = None
    offset = None

    # Now we can see if we have one good success:
    if len(success_dict) > 0:

        # We need to collapse entries with the same barcode here.
        # That is, it's possible that for different padded positions we're getting the same barcode back.
        # If this is the case, we should not call the correction ambiguous.

        lowest_dist = float("inf")
        collapsed_success_dict = dict()
        for lev_dist, match_info in success_dict.items():
            new_match_info = dict()

            # Arbitrarily keep only the first match:
            for barcode, status, index in match_info:
                if barcode not in new_match_info:
                    new_match_info[barcode] = (status, index)

            collapsed_success_dict[lev_dist] = [(barcode, status, index) for barcode, (status, index) in
                                                new_match_info.items()]

            # Since we're in here, let's just get the minimum distance now:
            if lev_dist < lowest_dist:
                lowest_dist = lev_dist

        # Now we can check that the lowest distance entry is a singleton:
        if len(collapsed_success_dict[lowest_dist]) == 1:
            new_bc = collapsed_success_dict[lowest_dist][0][0]
            result_status = collapsed_success_dict[lowest_dist][0][1]
            offset = collapsed_success_dict[lowest_dist][0][2]
        else:
            result_status = barcode_utils.SymSpellMatchResultType.AMBIGUOUS
    else:
        # Get the most prevalent failure type and report it:
        most_failures = -1
        for failure_type in failure_counts.keys():
            if failure_counts[failure_type] > most_failures:
                result_status = failure_type
                most_failures = failure_counts[failure_type]

    return new_bc, result_status, offset


def _get_barcode_tag_length_from_model(lb_model, barcode_tag):

    barcode_seg_name = lb_model.get_segment_name_for_annotation_tag(barcode_tag)
    if not barcode_seg_name:
        print(f"ERROR: Could not determine {lb_model.name} model segment from tag name: {barcode_tag}", file=sys.stderr)
        sys.exit(1)

    barcode_length = lb_model.get_segment_length(barcode_seg_name)
    if not barcode_length:
        print(f"ERROR: Could not extract {lb_model.name} model {barcode_seg_name} barcode segment from model.",
              file=sys.stderr)
        sys.exit(1)

    return barcode_length


def _generate_barcode_freqs_from_allow_list(barcode_allow_list):
    with tempfile.NamedTemporaryFile('w', prefix="tmp_barcode_allow_list_uniform_freqs",
                                     suffix=".tsv", dir=os.getcwd(), delete=False) as f:
        fname = f.name
        for bc in barcode_allow_list:
            f.write(f"{bc}\t1\n")

    return fname


def _filter_barcode_freqs_by_allow_list(barcode_freqs, barcode_allow_list):
    with open(barcode_freqs, 'r') as f:
        with tempfile.NamedTemporaryFile('w', prefix="tmp_barcode_freqs_filtered_to_allow_list",
                                         suffix=".tsv", dir=os.getcwd(), delete=False) as out_file:
            fname = out_file.name
            for line in f:
                bc = line.strip().split("\t")[0]
                if bc in barcode_allow_list:
                    out_file.write(line)
    return fname
