import logging
import sys
import itertools
import time
import os

import click
import tqdm

import pysam
import multiprocessing as mp

from ..utils import bam_utils
from ..utils import cli_utils


logger = logging.getLogger(__name__)


@click.command()
@cli_utils.input_pbi
@click.option(
    "-o",
    "--out-base-name",
    default="longbow_demultiplexed",
    required=False,
    show_default=True,
    type=str,
    help="base name for output files",
)
@click.option(
    "-d",
    "--demux-on-tag",
    type=str,
    default="YN",
    show_default=True,
    required=False,
    help="BAM tag on which to demultiplex (e.g. YN, BC)"
)
@cli_utils.input_bam
@click.pass_context
def main(ctx, pbi, out_base_name, demux_on_tag, input_bam):
    """Separate reads into files based on which model they fit best.

    Resulting reads will be annotated with the model they best fit as well as the score and segments for that model."""

    t_start = time.time()

    threads = ctx.obj["THREADS"]

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("Processing %d reads", read_count)
    else:
        read_count = bam_utils.get_read_count_from_bam_index(input_bam)
        if read_count:
            logger.info("Processing %d reads", read_count)

    # Create queues for data:
    queue_size = threads * 2 if threads < 10 else 20
    manager = mp.Manager()
    input_data_queue = manager.Queue(maxsize=queue_size)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_pool = []

    for i in range(threads):
        p = mp.Process(
            target=_worker_demux_fn, args=(input_data_queue, results, i, demux_on_tag)
        )
        p.start()
        worker_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file:

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header)

        # Start output worker:
        res = manager.dict({"num_reads_demultiplexed": 0, "num_reads_discarded": 0})
        output_worker = mp.Process(
            target=_write_thread_fn, args=(results, out_header, out_base_name, res, read_count)
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
        f"Demultiplexed {res['num_reads_demultiplexed']} reads. Discarded {res['num_reads_discarded']}."
    )

    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {(res['num_reads_demultiplexed']+res['num_reads_discarded'])/(et - t_start):2.2f} reads/s.")


def _write_thread_fn(out_queue, out_bam_header, out_bam_base_name, res, read_count):
    """Thread / process fn to write out all our data."""
    out_bam_header = pysam.AlignmentHeader.from_dict(out_bam_header)

    try:
        out_file_dict = {}

        with tqdm.tqdm(
            desc="Progress",
            unit=" read",
            colour="green",
            file=sys.stderr,
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
                read, tag_value = raw_data

                if tag_value is not None:
                    read = pysam.AlignedSegment.fromstring(read, out_bam_header)

                    if tag_value not in out_file_dict:
                        out_file_dict[tag_value] = pysam.AlignmentFile(
                            f'{out_bam_base_name}.{tag_value}.bam',
                            "wb",
                            header=out_bam_header
                        )

                    out_file_dict[tag_value].write(read)

                    res["num_reads_demultiplexed"] += 1
                else:
                    res["num_reads_discarded"] += 1

                pbar.update(1)
    finally:
        for out_file in out_file_dict.values():
            out_file.close()


def _worker_demux_fn(in_queue, out_queue, worker_num, demux_on_tag):
    """Function to run in each subthread / subprocess.
    Simply pulls user-specified tag from read."""

    num_reads_processed = 0

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

        tag_value = read.get_tag(demux_on_tag) if read.has_tag(demux_on_tag) else None

        out_queue.put((read.to_string(), tag_value))

        num_reads_processed += 1

    logger.debug(f"Worker {worker_num}: Num reads processed: {num_reads_processed}")
