import logging
import sys
import itertools
import re
import time
import os

import click
import click_log
import tqdm

import ssw

import pysam
import multiprocessing as mp

import gzip
from construct import *

from ..utils import bam_utils
from ..utils.bam_utils import SegmentInfo
from ..utils.model import LibraryModel


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("annotate")
click_log.basic_config(logger)


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
    type=str,
    default="mas15",
    show_default=True,
    help="The model to use for annotation.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, threads, output_bam, model, input_bam):
    """Annotate reads in a BAM file with segments from the model."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Get our model:
    if LibraryModel.has_prebuilt_model(model):
        logger.info(f"Using %s", LibraryModel.pre_configured_models[model]["description"])
        m = LibraryModel.build_pre_configured_model(model)
    else:
        logger.info(f"Loading model from json file: %s", model)
        m = LibraryModel.from_json_file(model)

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("Annotating %d reads", read_count)

    # Create queues for data:
    queue_size = threads * 2 if threads < 10 else 20
    manager = mp.Manager()
    input_data_queue = manager.Queue(maxsize=queue_size)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_pool = []

    for i in range(threads):
        p = mp.Process(
            target=_worker_segmentation_fn, args=(input_data_queue, results, i, m)
        )
        p.start()
        worker_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file:

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[m])

        # Start output worker:
        res = manager.dict({"num_reads_annotated": 0, "num_sections": 0})
        output_worker = mp.Process(
            target=_write_thread_fn,
            args=(results, out_header, output_bam, not sys.stdin.isatty(), res, read_count, m)
        )
        output_worker.start()

        # Add in a sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        iter_data = itertools.chain(bam_file, (None,) * threads)
        for r in iter_data:
            # We have to adjust for our sentinel value if we've got to it:
            if r is not None:
                r = r.to_string()
            input_data_queue.put(r)

        logger.debug("Finished reading data and sending it to sub-processes.")
        logger.debug("Waiting for sub-processes to finish...")

        # Wait for our input jobs to finish:
        for p in worker_pool:
            p.join()

        logger.debug("All workers stopped.")
        logger.debug("Terminating output process.")

        # Now that our input processes are done, we can add our exit sentinel onto the output queue and
        # wait for that process to end:
        results.put(None)
        output_worker.join()

    logger.info(
        f"Annotated {res['num_reads_annotated']} reads with {res['num_sections']} total sections."
    )
    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {res['num_reads_annotated']/(et - t_start):2.2f} reads/s.")


def get_segments(read):
    """Get the segments corresponding to a particular read by reading the segments tag information."""
    return read.to_string(), [
        SegmentInfo.from_tag(s) for s in read.get_tag(bam_utils.SEGMENTS_TAG).split(bam_utils.SEGMENT_TAG_DELIMITER)
    ]


def _write_thread_fn(out_queue, out_bam_header, out_bam_file_name, disable_pbar, res, read_count, model):
    """Thread / process fn to write out all our data."""

    with pysam.AlignmentFile(
        out_bam_file_name, "wb", header=out_bam_header
    ) as out_bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        colour="green",
        file=sys.stderr,
        disable=disable_pbar,
        total=read_count
    ) as pbar:

        ssw_aligner = ssw.Aligner()

        while True:
            # Wait for some output data:
            raw_data = out_queue.get()

            # Check for exit sentinel:
            if raw_data is None:
                break
            # Should really never be None, but just in case:
            elif raw_data is None:
                continue

            # Unpack data:
            read, ppath, logp, is_rc = raw_data

            # Condense the output annotations so we can write them out with indices:
            segments = bam_utils.collapse_annotations(ppath)

            read = pysam.AlignedSegment.fromstring(read, out_bam_header)

            # Obligatory log message:
            logger.debug(
                "Path for read %s (%2.2f)%s: %s",
                read.query_name,
                logp,
                " (RC)" if is_rc else "",
                segments,
            )

            # Write our our read:
            bam_utils.write_annotated_read(read, segments, is_rc, logp, model, ssw_aligner, out_bam_file)

            # Increment our counters:
            res["num_reads_annotated"] += 1
            res["num_sections"] += len(segments)

            pbar.update(1)


def _worker_segmentation_fn(in_queue, out_queue, worker_num, model):
    """Function to run in each subthread / subprocess.
    Segments each read and place the segments in the output queue."""

    num_reads_segmented = 0

    while True:
        # Wait until we get some data.
        # Note: Because we have a sentinel value None inserted at the end of the input data for each
        #       subprocess, we don't have to add a timeout - we're guaranteed each process will always have
        #       at least one element.
        raw_data = in_queue.get()

        # Check for exit sentinel:
        if raw_data is None:
            break
        # Should really never be None, but just in case:
        elif raw_data is None:
            continue

        # Unpack our data here:
        read = raw_data
        read = pysam.AlignedSegment.fromstring(
            read, pysam.AlignmentHeader.from_dict(dict())
        )

        # Process and place our data on the output queue:
        segment_info = _segment_read(read, model)

        out_queue.put(segment_info)
        num_reads_segmented += 1

    logger.debug(f"Worker %d: Num reads segmented: %d", worker_num, num_reads_segmented)


def _segment_read(read, model):
    is_rc = False
    logp, ppath = model.annotate(read.query_sequence)

    rc_logp, rc_ppath = model.annotate(bam_utils.reverse_complement(read.query_sequence))
    if rc_logp > logp:
        logp = rc_logp
        ppath = rc_ppath
        is_rc = True
        logger.debug("Sequence scored better in RC: %s", read.query_name)

    return read.to_string(), ppath, logp, is_rc
