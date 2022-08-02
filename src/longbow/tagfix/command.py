import logging
import re
import time
import sys
import itertools

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp

from construct import *

import longbow.utils.constants
from ..utils import bam_utils
from ..utils.cli_utils import format_obnoxious_warning_message


PROG_NAME = "tagfix"

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(PROG_NAME)
click_log.basic_config(logger)


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
    "--output-bam",
    default="-",
    type=click.Path(exists=False),
    help="annotated bam output  [default: stdout]",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(threads, output_bam, force, input_bam):
    """Update longbow read tags after alignment."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Configure process manager:
    # NOTE: We're using processes to overcome the Global Interpreter Lock.
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(threads)
    results = manager.Queue()

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

        # Start worker sub-processes:
        res = manager.dict({"num_tags_corrected": 0, "num_reads": 0})
        worker_process_pool = []
        for _ in range(threads):
            p = mp.Process(
                target=_correct_read_tags, args=(process_input_data_queue, results, bam_file.header, res)
            )
            p.start()
            worker_process_pool.append(p)

        # Check if we have run this command on a bam file already:
        if bam_utils.bam_header_has_longbow_command_program_group(bam_file.header, PROG_NAME):
            logger.error(f"{PROG_NAME} has already been run on input bam: {input_bam}")
            sys.exit(1)

        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header)

        # Start output worker:
        output_worker = mp.Process(
            target=_output_writer_fn,
            args=(
                results,
                out_header,
                output_bam,
                pbar,
            ),
        )
        output_worker.start()

        # Add in a sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        iter_data = itertools.chain(bam_file, (None,) * threads)
        for r in iter_data:
            if r is not None:
                process_input_data_queue.put(r.to_dict())
            else:
                process_input_data_queue.put(r)

        # Wait for our input jobs to finish:
        for p in worker_process_pool:
            p.join()

        # Now that our input processes are done, we can add our exit sentinel onto the output queue and
        # wait for that process to end:
        results.put(None)
        output_worker.join()

    logger.info(f"Corrected tags in {res['num_tags_corrected']} reads of {res['num_reads']} total.")
    
    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {res['num_reads']/(et - t_start):2.2f} reads/s.")

    if res['num_tags_corrected'] == 0:
        logger.warning(format_obnoxious_warning_message("No read tags were corrected.  This is very unlikely.  "
                                                        "You should check your data."))


def _output_writer_fn(out_queue, out_bam_header, out_bam_file_name, pbar):
    """Thread / process fn to write out all our data."""

    with pysam.AlignmentFile(
        out_bam_file_name, "wb", header=out_bam_header
    ) as out_bam_file:
        while True:
            # Wait for some output data:
            raw_data = out_queue.get()

            # Check for exit sentinel:
            if raw_data is None:
                break

            # Unpack data:
            read = raw_data
            read = pysam.AlignedSegment.from_dict(read, out_bam_header)

            # Write out our read:
            out_bam_file.write(read)

            pbar.update(1)


def _correct_read_tags(in_queue, out_queue, bam_header, res):
    """Function to run in each subprocess.
    Do the tag correction here."""

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
        read = pysam.AlignedSegment.from_dict(raw_data, bam_header)

        if read.is_reverse:

            # Reverse the segments and update the positions:
            segments_tag = read.get_tag(longbow.utils.constants.SEGMENTS_TAG)
            segments = segments_tag.split(longbow.utils.constants.SEGMENT_TAG_DELIMITER)
            read_length = len(read.query_sequence)
            new_segments = []
            for seg in reversed(segments):
                seg_name, start, end = re.split("[:-]", seg)
                new_start = read_length - int(end)
                new_end = read_length - int(start)

                new_segments.append(f"{seg_name}:{new_start}-{new_end}")

            read.set_tag(longbow.utils.constants.SEGMENTS_TAG,
                         longbow.utils.constants.SEGMENT_TAG_DELIMITER.join(new_segments))

            # Reverse the Segment Qualities:
            seg_quals = read.get_tag(longbow.utils.constants.SEGMENTS_QUAL_TAG)
            new_seg_quals = seg_quals.split(longbow.utils.constants.SEGMENT_TAG_DELIMITER)[::-1]
            read.set_tag(longbow.utils.constants.SEGMENTS_QUAL_TAG,
                         longbow.utils.constants.SEGMENT_TAG_DELIMITER.join(new_seg_quals))

            # Increment our counter:
            res["num_tags_corrected"] += 1

        out_queue.put(read.to_dict())
        res["num_reads"] += 1
