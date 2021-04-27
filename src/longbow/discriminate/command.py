import logging
import sys
import itertools
import re
import time
import collections
import os
import math

from inspect import getframeinfo, currentframe, getdoc

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp

from ..utils import bam_utils
from ..utils.model import reverse_complement
from ..utils.model import LibraryModel

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("discriminate")
click_log.basic_config(logger)


# Named tuple to store alignment information:
class SegmentInfo(collections.namedtuple("SegmentInfo", ["name", "start", "end"])):

    _tag_regex = re.compile(r"(.*?):(\d+)-(\d+)")

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return f"SegmentInfo({self.to_tag()})"

    def to_tag(self):
        return f"{self.name}:{self.start}-{self.end}"

    @classmethod
    def from_tag(cls, tag_string):
        match = cls._tag_regex.match(tag_string)
        return SegmentInfo(match[1], int(match[2]), int(match[3]))


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
    "--out-base-name",
    default="longbow_discriminated",
    required=False,
    show_default=True,
    type=str,
    help="base name for output files",
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, out_base_name, threads, input_bam):
    """Separate reads into files based on which model they fit best.

    Resulting reads will be annotated with the model they best fit as well as the score and segments for that model."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

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

    # Create a dictionary of models to use to annotate and score our reads:
    model_dict = {
        "mas15": LibraryModel.build_and_return_mas_seq_model(),
        "mas10": LibraryModel.build_and_return_mas_seq_10_model(),
    }

    # Start worker sub-processes:
    worker_pool = []

    for i in range(threads):
        p = mp.Process(
            target=_worker_segmentation_fn, args=(input_data_queue, results, model_dict, i)
        )
        p.start()
        worker_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file:

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group("discriminate", bam_file.header)

        # Start output worker:
        res = manager.dict({"num_reads_annotated": 0, "num_sections": 0,
                            "model_reads_annotated_count_dict": dict()})
        output_worker = mp.Process(
            target=_write_thread_fn, args=(results, out_header, out_base_name, model_dict, res, read_count)
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
    for model_name, count in res["model_reads_annotated_count_dict"].items():
        logger.info(f"Model {model_name} annotated {count} reads.")

    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {res['num_reads_annotated']/(et - t_start):2.2f} reads/s.")


def _write_thread_fn(out_queue, out_bam_header, out_bam_base_name, model_dict, res, read_count):
    """Thread / process fn to write out all our data."""

    try:
        out_file_dict = {name: pysam.AlignmentFile(f"{out_bam_base_name}_{name}.bam", "wb", header=out_bam_header) for
                         name in model_dict.keys()}
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
                read, ppath, logp, is_rc, model_name = raw_data
                read = pysam.AlignedSegment.fromstring(read, out_bam_header)

                # Condense the output annotations so we can write them out with indices:
                segments = _collapse_annotations(ppath)

                # Obligatory log message:
                logger.debug(
                    "Path for read %s (%2.2f)%s: %s",
                    read.query_name,
                    logp,
                    " (RC)" if is_rc else "",
                    segments,
                )

                # Set our tag and write out the read to the annotated file:
                read.set_tag(
                    bam_utils.SEGMENTS_TAG, bam_utils.SEGMENT_TAG_DELIMITER.join([s.to_tag() for s in segments])
                )

                # If we're reverse complemented, we make it easy and just reverse complement the read and add a
                # tag saying that the read was RC:
                read.set_tag(bam_utils.SEGMENTS_RC_TAG, is_rc)
                read.set_tag(bam_utils.READ_MODEL_SCORE_TAG, logp)
                read.set_tag(bam_utils.READ_MODEL_NAME_TAG, model_name)
                if is_rc:
                    quals = read.query_qualities[::-1]
                    seq = reverse_complement(read.query_sequence)
                    read.query_sequence = seq
                    read.query_qualities = quals
                out_file_dict[model_name].write(read)

                # Increment our counters:
                res["num_reads_annotated"] += 1
                res["num_sections"] += len(segments)
                count_dict = res["model_reads_annotated_count_dict"]
                try:
                    count_dict[model_name] += 1
                except KeyError:
                    count_dict[model_name] = 1
                res["model_reads_annotated_count_dict"] = count_dict

                pbar.update(1)
    finally:
        for out_file in out_file_dict.values():
            out_file.close()


def _worker_segmentation_fn(in_queue, out_queue, model_dict, worker_num):
    """Function to run in each subthread / subprocess.
    Annotates each read with both models and assigns the read to the model with the better score."""

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
        segment_info = _annotate_and_assign_read_to_model(read, model_dict)

        out_queue.put(segment_info)
        num_reads_segmented += 1

    logger.debug(f"Worker %d: Num reads segmented: %d", worker_num, num_reads_segmented)


def _collapse_annotations(path):
    """Collapses given path into a list of SegmentInfo objects."""
    last = ""
    start = 0
    segments = []
    i = 0
    for i, seg in enumerate(path):
        if seg != last:
            if i != 0:
                segments.append(SegmentInfo(last, start, i - 1))
            last = seg
            start = i
    # Don't forget the last one:
    segments.append(SegmentInfo(last, start, i))

    return segments


def _annotate_and_assign_read_to_model(read, model_dict):
    """Annotate the given read with all given models and assign the read to the model with the best score."""

    best_model = ""
    best_logp = -math.inf
    best_path = None
    best_fit_is_rc = False
    model_scores = dict()
    for name, model in model_dict.items():

        _, ppath, logp, is_rc = _annotate_read(read, model)

        model_scores[name] = logp

        if logp > best_logp:
            best_model = name
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

    rc_logp, rc_ppath = model.annotate(reverse_complement(read.query_sequence))
    if rc_logp > logp:
        logp = rc_logp
        ppath = rc_ppath
        is_rc = True

    return read.to_string(), ppath, logp, is_rc
