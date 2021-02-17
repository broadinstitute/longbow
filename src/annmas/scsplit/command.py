import logging
import sys
import itertools

import time

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp

from inspect import getframeinfo, currentframe, getdoc

from ..utils.model import array_element_structure

from ..annotate.command import SegmentInfo
from ..annotate.command import SEGMENTS_TAG
from ..annotate.command import _get_segments

from ..utils.constants import __SENTINEL_VALUE

from ..meta import VERSION

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("scsplit")
click_log.basic_config(logger)


__DEFAULT_DUMMY_CBC = "CTGCCTAACCTGATCC"
__DEFAULT_OUT_BASE_NAME = "scsplit"
__DEFAULT_UMI_LENGTH = 10


__OUT_READ_FILE_SUFFIX = "_mates"
__OUT_WHITELIST_FILE_SUFFIX = "_whitelist.txt"


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
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
    "--output-base-name",
    default=__DEFAULT_OUT_BASE_NAME,
    type=str,
    help=f"base name for output files [default: {__DEFAULT_OUT_BASE_NAME}]",
)
@click.option(
    "-c",
    "--cell-barcode",
    default=__DEFAULT_DUMMY_CBC,
    type=str,
    help=f"dummy cell barcode to use for the dataset [default: {__DEFAULT_DUMMY_CBC}, "
         f"length: {len(__DEFAULT_DUMMY_CBC)}]",
)
@click.option(
    "-u",
    "--umi-length",
    default=__DEFAULT_UMI_LENGTH,
    type=int,
    show_default=True,
    help=f"length of the UMI from this library",
)
@click.argument(
    "input-bam",
    default="-" if not sys.stdin.isatty() else None,
    type=click.File("rb"),
)
def main(threads, output_base_name, cell_barcode, umi_length, input_bam):
    f"""Create a set of files suitable for use in `alevin` for single-cell analysis.
    This tool coerces a set of reads from a single source into a format that `alevin` can ingest.
    
    Segment names are assumed to be those in the default model (utils/model.py).
    
    INPUT_BAM should contain reads that have been processed by `annmas segment`.

    The output from this tool consists of several files:
        {__DEFAULT_OUT_BASE_NAME}{__OUT_READ_FILE_SUFFIX}1.fastq: 
            A file containing partial sequences for all reads in the given input file.  These partial reads consist of the 
            dummy cell barcode + detected UMI for each read in the given input file.
        {__DEFAULT_OUT_BASE_NAME}{__OUT_READ_FILE_SUFFIX}2.fastq: 
            A file containing partial sequences for all reads in the given input file.  These partial reads consist of the 
            transcript sequences for all reads in the given input file.  Transcript sequences include data after the UMI 
            and before the Poly-A tail.  All bases outside of this range are excluded from the output.
        {__DEFAULT_OUT_BASE_NAME}{__OUT_WHITELIST_FILE_SUFFIX}: 
            A whitelist file for alevin containing the given dummy cell barcode. 
    """

    t_start = time.time()

    logger.info("Invoked via: annmas %s", " ".join(sys.argv[1:]))

    # Make sure we're given an input bam file we can work with:
    _validate_input_bam(input_bam)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Configure process manager:
    # NOTE: We're using processes to overcome the Global Interpreter Lock.
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(threads)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_process_pool = []
    for _ in range(threads):
        p = mp.Process(
            target=_sub_process_work_fn, args=(process_input_data_queue, results, umi_length)
        )
        p.start()
        worker_process_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        colour="green",
        file=sys.stderr,
        leave=False,
        disable=not sys.stdin.isatty(),
    ) as pbar:

        # Start output worker:
        res = manager.dict({"num_reads_processed": 0})
        output_worker = mp.Process(
            target=_sub_process_write_fn,
            args=(
                results,
                output_base_name,
                cell_barcode,
                pbar,
                res,
            ),
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

    # Write out our CBC whitelist file:
    with open(f"{output_base_name}{__OUT_WHITELIST_FILE_SUFFIX}", "w") as f:
        f.write(f"{cell_barcode}\n")

    logger.info(f"Processed {res['num_reads_segmented']} reads.")
    logger.info(f"CBC length: {len(cell_barcode)}.")
    logger.info(f"UMI length: {umi_length}.")
    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)


def _validate_input_bam(input_bam_path):
    """Assert that the given input_bam was created by `annmas segment`."""
    with pysam.AlignmentFile(
            input_bam_path, "rb", check_sq=False, require_index=False
    ) as bam_file:
        in_bam_header_dict = bam_file.header.to_dict()
        if "PG" not in in_bam_header_dict:
            logger.warn("Could not find PG entry in header.  Cannot confirm that this file is compatible.")
        else:
            found_segment_cmd = False
            for info in [item for item in in_bam_header_dict["PG"] if item["PN"] == "annmas"]:
                if info["ID"].split("-")[1] == "segment":
                    found_segment_cmd = True
                    break
            if not found_segment_cmd:
                logger.error(
                    "Input bam file header does not indicate that it was created by annmas segment.  This tool requires `annmas segment` reads as input data.")
                sys.exit(1)


def _get_named_segment_from_list(seg_list, name, read_name):
    """Get the segment with the given name from the list of SegmentInfo objects.
    If multiple items have the given name a warning is logged and None is returned."""

    found_segments = [s for s in seg_list if s.name == name]
    if len(found_segments) != 1:
        logger.warn("Could not process read: %s - multiple %s annotations are present (%d)",
                    read_name, name, len(found_segments))
        return None
    return found_segments[0]


def _sub_process_work_fn(in_queue, out_queue, umi_length):
    """Function to run in each subprocess.
    Extracts and returns all segments from an input read."""
    while True:
        # Wait until we get some data.
        # Note: Because we have a sentinel value None inserted at the end of the input data for each
        #       subprocess, we don't have to add a timeout - we're guaranteed each process will always have
        #       at least one element.
        raw_data = in_queue.get()

        # Check for exit sentinel:
        if raw_data is None:
            return

        # Unpack our data here:
        read = pysam.AlignedSegment.fromstring(
            raw_data, pysam.AlignmentHeader.from_dict(dict())
        )
        _, segments = _get_segments(read)

        # Get the 10x adapter position:
        tenx_adapter_segment = _get_named_segment_from_list(segments, "10x_Adapter", read.query_name)
        if tenx_adapter_segment is None:
            continue

        # Get the Poly-A position:
        poly_a_segment = _get_named_segment_from_list(segments, "Poly_A", read.query_name)
        if poly_a_segment is None:
            continue

        # Now we grab the bases just after the 10x adapter as the UMI
        # and the bases between the UMI and the poly A for the transcript

        # Note: Positions are inclusive so we must add 1 to the end position to get that base as well:
        umi_start = tenx_adapter_segment.end+1
        umi_end = umi_start + umi_length
        umi_bases = read.query_sequence[umi_start:umi_end+1]
        umi_quals = read.query_alignment_qualities[umi_start:umi_end+1]

        transcript_bases = read.query_sequence[umi_end+1:poly_a_segment.start]
        transcript_quals = read.query_alignment_qualities[umi_end+1:poly_a_segment.start]

        # Place our data on the output queue:
        out_queue.put(tuple(read.query_name, umi_bases, umi_quals, transcript_bases, transcript_quals))


def _sub_process_write_fn(
    out_queue,
    out_base_name,
    cell_barcode,
    pbar,
    res,
):
    """Thread / process fn to write out all our data."""

    with open(f"{out_base_name}{__OUT_READ_FILE_SUFFIX}1.fastq", "w") as mates1_file, \
         open(f"{out_base_name}{__OUT_READ_FILE_SUFFIX}2.fastq", "w") as mates2_file:

        while True:
            # Wait for some output data:
            raw_data = out_queue.get()

            # Check for exit sentinel:
            if raw_data is None:
                break

            # Unpack data:
            read_name, umi_bases, umi_quals, transcript_bases, transcript_quals = raw_data

            # Create mates1 and mates2 records:
            mates_1_record = pysam.FastxRecord(
                name=read_name,
                sequence=cell_barcode + umi_bases,
                quality=(chr(33 + 60) * len(cell_barcode)) + umi_quals
            )
            mates_2_record = pysam.FastxRecord(
                name=read_name,
                sequence=transcript_bases,
                quality=transcript_quals
            )

            # Write out mates1 and mates2 records:
            mates1_file.write(str(mates_1_record))
            mates2_file.write(str(mates_2_record))

            # Increment our counters:
            res["num_reads_processed"] += 1
            pbar.update(1)

            # Obligatory log message:
            logger.debug("Processed read: %s", read_name)
