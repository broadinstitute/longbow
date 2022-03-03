import logging
import math
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

import longbow.utils.constants
from ..utils import bam_utils
from ..utils.bam_utils import SegmentInfo
from ..utils import model as LongbowModel
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
    multiple=True,
    show_default=True,
    help="The model(s) to use for annotation.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
)
@click.option(
    "-c",
    "--chunk",
    type=str,
    default="",
    required=False,
    help="Process a single chunk of data (e.g. specify '2/4' to process the second of four equally-sized "
         "chunks across the dataset)"
)
@click.option(
    "-l",
    "--min-length",
    type=int,
    default=0,
    show_default=True,
    required=False,
    help="Minimum length of a read to process.  Reads shorter than this length will not be annotated."
)
@click.option(
    "-L",
    "--max-length",
    type=int,
    default=longbow.utils.constants.DEFAULT_MAX_READ_LENGTH,
    show_default=True,
    required=False,
    help="Maximum length of a read to process.  Reads longer than this length will not be annotated."
)
@click.option(
    "--min-rq",
    type=float,
    default=-2,
    show_default=True,
    required=False,
    help="Minimum ccs-determined read quality for a read to be annotated.  CCS read quality range is [-1,1]."
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
def main(pbi, threads, output_bam, model, chunk, min_length, max_length, min_rq, force, input_bam):
    """Annotate reads in a BAM file with segments from the model."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Get our model(s):
    lb_models = bam_utils.load_models(model, input_bam)

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    read_num = 0
    start_offset = 0
    end_offset = math.inf

    if not os.path.exists(pbi) and chunk is not "":
        raise ValueError(f"Chunking specified but pbi file '{pbi}' not found")

    if os.path.exists(pbi):
        if chunk is not "":
            (chunk, num_chunks) = re.split("/", chunk)
            chunk = int(chunk)
            num_chunks = int(num_chunks)

            # Decode PacBio .pbi file and determine the shard offsets.
            offsets, zmw_counts, read_count, read_counts_per_chunk, read_nums = bam_utils.compute_shard_offsets(pbi, num_chunks)

            start_offset = offsets[chunk - 1]
            end_offset = offsets[chunk] if chunk < len(offsets) else offsets[chunk - 1]
            read_count = read_counts_per_chunk[chunk - 1] if chunk < len(offsets) else 0
            read_num = read_nums[chunk - 1] if chunk < len(offsets) else 0

            logger.info("Annotating %d reads from chunk %d/%d (reads %d-%d)", read_count, chunk, num_chunks, read_num, read_num + read_count - 1)
        else:
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
            target=_worker_segmentation_fn,
            args=(input_data_queue, results, i, lb_models, min_length, max_length, min_rq)
        )
        p.start()
        worker_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam if start_offset == 0 else input_bam.name, "rb", check_sq=False, require_index=False
    ) as bam_file:

        # If we're chunking, advance to the specified virtual file offset.
        if start_offset > 0:
            bam_file.seek(start_offset)

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=lb_models)

        # Start output worker:
        res = manager.dict({"num_reads_annotated": 0, "num_sections": 0})
        output_worker = mp.Process(
            target=_write_thread_fn,
            args=(results, out_header, output_bam, not sys.stdin.isatty(), res, read_count, lb_models)
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

            if start_offset > 0:
                if bam_file.tell() >= end_offset or r is None:
                    [input_data_queue.put(None) for _ in range(threads)]
                    break

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
        SegmentInfo.from_tag(s) for s in read.get_tag(longbow.utils.constants.SEGMENTS_TAG).split(
            longbow.utils.constants.SEGMENT_TAG_DELIMITER)
    ]


def _write_thread_fn(out_queue, out_bam_header, out_bam_file_name, disable_pbar, res, read_count, lb_models):
    """Thread / process fn to write out all our data."""

    lb_models_dict = {k.name: k for k in lb_models}

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
            read, ppath, logp, is_rc, model_name = raw_data

            if read is not None:
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
                bam_utils.write_annotated_read(read, segments, is_rc, logp, lb_models_dict[model_name], ssw_aligner, out_bam_file)

                # Increment our counters:
                res["num_reads_annotated"] += 1
                res["num_sections"] += len(segments)

            pbar.update(1)


def _worker_segmentation_fn(in_queue, out_queue, worker_num, lb_models, min_length, max_length, min_rq):
    """Function to run in each subthread / subprocess.
    Segments each read and place the segments in the output queue."""

    num_reads_processed, num_reads_segmented = 0, 0

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

        # Check for min/max length and min quality:
        segment_info = None, None, None, None, None
        if len(read.query_sequence) < min_length:
            logger.debug(f"Read is shorter than min length.  "
                         f"Skipping: {read.query_name} ({len(read.query_sequence)} < {min_length})")
        elif len(read.query_sequence) > max_length:
            logger.debug(f"Read is longer than max length.  "
                         f"Skipping: {read.query_name} ({len(read.query_sequence)} > {max_length})")
        elif read.get_tag("rq") < min_rq:
            logger.debug(f"Read quality is below the minimum.  "
                         f"Skipping: {read.query_name} ({read.get_tag('rq')} < {min_rq})")
        else:
            # Process and place our data on the output queue:
            segment_info = _annotate_and_assign_read_to_model(read, lb_models)
            num_reads_segmented += 1

        out_queue.put(segment_info)
        num_reads_processed += 1

    logger.debug(f"Worker %d: Num reads segmented/processed: %d/%d", worker_num, num_reads_segmented, num_reads_processed)


def _annotate_and_assign_read_to_model(read, model_list):
    """Annotate the given read with all given models and assign the read to the model with the best score."""

    best_model = ""
    best_logp = -math.inf
    best_path = None
    best_fit_is_rc = False
    model_scores = dict()
    for model in model_list:

        _, ppath, logp, is_rc = _annotate_read(read, model)

        model_scores[model.name] = logp

        if logp > best_logp:
            best_model = model.name
            best_logp = logp
            best_path = ppath
            best_fit_is_rc = is_rc

    # Provide some info as to which model was chosen:
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug("%s model scores: %s", read.query_name, str(model_scores))
        logger.debug(
            "Sequence %s scored best with model%s: %s (%2.4f)",
            read.query_name,
            " in RC " if best_fit_is_rc else "",
            best_model,
            best_logp
        )

    return read.to_string(), best_path, best_logp, best_fit_is_rc, best_model


def _annotate_read(read, model):
    is_rc = False
    logp, ppath = model.annotate(read.query_sequence)

    rc_logp, rc_ppath = model.annotate(bam_utils.reverse_complement(read.query_sequence))
    if rc_logp > logp:
        logp = rc_logp
        ppath = rc_ppath
        is_rc = True

    return read.to_string(), ppath, logp, is_rc