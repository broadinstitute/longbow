from cmath import inf
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

import numpy as np
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
logger = logging.getLogger("peek")
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
    "--output-model",
    default="-",
    type=click.Path(exists=False),
    help="model name output  [default: stdout]",
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
    "-n",
    "--num-reads",
    type=int,
    default=100,
    required=False,
    help="Number of reads to examine when trying to guess the array model"
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
@click.option(
    "-d",
    "--include-deprecated-models",
    is_flag=True,
    default=False,
    show_default=True,
    help="Examine the deprecated built-in models as well",
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, threads, output_model, chunk, num_reads, min_length, max_length, min_rq, force, include_deprecated_models, input_bam):
    """Guess the best pre-built array model to use for annotation."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_model, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Make all prebuilt models
    models = {}
    for model_name in LibraryModel.pre_configured_models:
        if not LibraryModel.pre_configured_models[model_name]['deprecated'] or include_deprecated_models:
            m = LibraryModel.build_pre_configured_model(model_name)
            models[model_name] = m

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

            logger.info("Detecting best model in %d reads from chunk %d/%d (reads %d-%d)", min(read_count, num_reads), chunk, num_chunks, read_num, read_num + read_count - 1)
        else:
            read_count = bam_utils.load_read_count(pbi)
            logger.info("Detecting best model in %d reads", min(read_count, num_reads))
    else:
        logger.info("Detecting best model in %d reads", num_reads)

    # Create queues for data:
    queue_size = threads * 2 if threads < 10 else 20
    manager = mp.Manager()
    input_data_queue = manager.Queue(maxsize=queue_size)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_pool = []

    for i in range(threads):
        p = mp.Process(
            target=_worker_attempt_segmentations_fn, args=(input_data_queue, results, i, min_length, max_length, min_rq, models)
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

        # Start output worker:
        res = manager.dict()
        output_worker = mp.Process(
            target=_collect_thread_fn,
            args=(results, output_model, not sys.stdin.isatty(), res, num_reads if read_count is None else min(read_count, num_reads))
        )
        output_worker.start()

        # Add in a sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        reads_seen = 0
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

            if reads_seen >= num_reads:
                [input_data_queue.put(None) for _ in range(threads)]
                break

            # Increment reads seen
            reads_seen += 1

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

    logger.info("Histogram of most likely model counts:")

    res = dict(sorted(res.items(), key=lambda item: item[1], reverse=True))
    plot_model_counts(res, models)
    best_model, best_model_count = next(iter(res.items()))

    logger.info(f"Overall most likely model: {best_model} (seen in {best_model_count} reads, {100.0*best_model_count/np.sum(list(res.values())):.1f}%)")

    with open(output_model if output_model != "-" else "/dev/stdout", "w") as wm:
        wm.write(f'{best_model}\n')

    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {reads_seen/(et - t_start):2.2f} reads/s.")


def plot_model_counts(res, models, max_width=50.0):
    big_block_char = u"\u2588"
    small_block_char = u"\u258F"

    model_tot = np.sum(list(res.values()))
    max_label_width = np.max(list(map(lambda x: len(x), models.keys())))
    for model_name in res:
        pct = 100.0 * res[model_name] / model_tot
        num_blocks = int(res[model_name] * max_width / model_tot)

        logger.info(f"  {model_name:>{max_label_width}} {big_block_char*(num_blocks+1)} {res[model_name]} ({pct:.1f}%)")

    for model_name in models:
        if model_name not in res:
            logger.info(f"  {model_name:>{max_label_width}} {small_block_char} 0 (0.0%)")


def get_segments(read):
    """Get the segments corresponding to a particular read by reading the segments tag information."""
    return read.to_string(), [
        SegmentInfo.from_tag(s) for s in read.get_tag(longbow.utils.constants.SEGMENTS_TAG).split(
            longbow.utils.constants.SEGMENT_TAG_DELIMITER)
    ]


def _collect_thread_fn(out_queue, out_bam_file_name, disable_pbar, res, read_count):
    """Thread / process fn to write out all our data."""

    with tqdm.tqdm(
            desc="Progress",
            unit=" read",
            colour="green",
            file=sys.stderr,
            disable=disable_pbar,
            total=read_count
        ) as pbar:

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
            best_model = raw_data

            if best_model not in res:
                res[best_model] = 0
            res[best_model] += 1

            pbar.update(1)


def _worker_attempt_segmentations_fn(in_queue, out_queue, worker_num, min_length, max_length, min_rq, models):
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

        # Check for min/max length and min quality:
        if len(read.query_sequence) < min_length:
            logger.debug(f"Read is shorter than min length.  "
                           f"Skipping: {read.query_name} ({len(read.query_sequence)} < {min_length})")
            continue
        elif len(read.query_sequence) > max_length:
            logger.debug(f"Read is longer than max length.  "
                           f"Skipping: {read.query_name} ({len(read.query_sequence)} > {max_length})")
            continue
        elif read.get_tag("rq") < min_rq:
            logger.debug(f"Read quality is below the minimum.  "
                           f"Skipping: {read.query_name} ({read.get_tag('rq')} < {min_rq})")
            continue

        # Process and place our data on the output queue:
        best_model = _get_best_model(read, models)

        out_queue.put(best_model)
        num_reads_segmented += 1

    logger.debug(f"Worker %d: Num reads segmented: %d", worker_num, num_reads_segmented)


def _get_best_model(read, models):
    model_scores = {}

    for model_name in models:
        model = models[model_name]
        fw_logp, _ = model.annotate(read.query_sequence)
        rc_logp, _ = model.annotate(bam_utils.reverse_complement(read.query_sequence))

        model_scores[f'{model_name}:fw'] = fw_logp
        model_scores[f'{model_name}:rc'] = rc_logp

    model_scores = dict(sorted(model_scores.items(), key=lambda item: item[1], reverse=True))

    best_model_full = next(iter(model_scores))
    best_model = re.sub(':.*$', '', best_model_full)

    return best_model
