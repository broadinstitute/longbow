import logging
import sys
import itertools
import math
import re

from collections import namedtuple

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp
import threading

from ..utils.model import build_default_model
from ..utils.model import reverse_complement
from ..utils.model import annotate

from ..meta import VERSION


# Named tuple to store alignment information:
class SegmentInfo(namedtuple("SegmentInfo", ["name", "start", "end"])):

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


SEGMENTS_TAG = "SG"
SEGMENTS_RC_TAG = "RC"

logger = logging.getLogger(__name__)
click_log.basic_config(logger)

mod_name = "annotate"


@click.command(name=mod_name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-m",
    "--model",
    required=False,
    type=click.Path(exists=True),
    help="pre-trained model to apply",
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
    required=True,
    type=click.Path(exists=False),
    help="segment-annotated bam output",
)
@click.argument("input-bam", type=click.Path(exists=True))
def main(model, threads, output_bam, input_bam):
    """Apply annotation and segmentation model to BAM file"""
    logger.info(f"annmas: {mod_name} started")

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} sub-processor(s)")

    m = build_default_model()
    if model is not None:
        m.from_yaml(model)
        logger.info(f"Using pretrained annotation model {model}")
    else:
        logger.info("Using default annotation model")

    # Configure process manager:
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(threads)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_process_pool = []
    for _ in range(threads):
        p = mp.Process(target=_sub_process_work_fn, args=(process_input_data_queue, results))
        p.start()
        worker_process_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file:

        with tqdm.tqdm(
            desc="Progress", unit=" read", colour="green", file=sys.stdout
        ) as pbar:

            # Get our header from the input bam file:
            out_bam_header_dict = bam_file.header.to_dict()

            # Add our program group to it:
            pg_dict = {
                "ID": f"annmas-annotate-{VERSION}",
                "PN": "annmas",
                "VN": f"{VERSION}",
                "DS": "Apply segmentation/annotation model to BAM file.",
                "CL": " ".join(sys.argv),
            }
            if "PG" in out_bam_header_dict:
                out_bam_header_dict["PG"].append(pg_dict)
            else:
                out_bam_header_dict["PG"] = [pg_dict]

            out_header = pysam.AlignmentHeader.from_dict(out_bam_header_dict)

            # Start output thread:
            # NOTE: we use a thread here so we can reuse the read objects without serializing them:
            # output_worker = threading.Thread(
            output_worker = mp.Process(
                target=_sub_thread_write_fn,
                args=(results, out_header, output_bam, pbar)
            )
            output_worker.start()

            # Add in a `None` sentinel value at the end of the queue - one for each subprocess - so we guarantee
            # that all subprocesses will exit:
            iter_data = itertools.chain(bam_file, (None,) * threads)
            for r in iter_data:
                if r is not None:
                    process_input_data_queue.put(
                        (r.to_string(), m)
                    )
                else:
                    process_input_data_queue.put(r)

            # Wait for our input jobs to finish:
            for p in worker_process_pool:
                p.join()

            # Now that our input processes are done, we can add our exit sentinel onto the output queue and
            # wait for that process to end:
            results.put(None)
            output_worker.join()

    logger.info(f"annmas: {mod_name} finished.")


def _sub_thread_write_fn(out_queue, out_bam_header, out_bam_file_name, pbar):
    """Thread / process fn to write out all our data."""

    num_reads_annotated = 0
    num_sections = 0

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
            read, ppath, logp, is_rc = raw_data
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
            read.set_tag(SEGMENTS_TAG, "|".join([s.to_tag() for s in segments]))

            # If we're reverse complemented, we make it easy and just reverse complement the read and add a tag saying
            # that the read was RC:
            read.set_tag(SEGMENTS_RC_TAG, is_rc)
            if is_rc:
                quals = read.query_qualities[::-1]
                seq = reverse_complement(read.query_sequence)
                read.query_sequence = seq
                read.query_qualities = quals
            out_bam_file.write(read)

            # Increment our counters:
            num_reads_annotated += 1
            num_sections += len(segments)

            pbar.update(1)

    logger.info(
        f"annmas {mod_name}: annotated {num_reads_annotated} reads with {num_sections} total sections."
    )
    return


def _sub_process_work_fn(in_queue, out_queue):
    """Function to run in each subprocess.
    Will segment each read and place the segments in the output queue."""
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
        read, um = raw_data
        read = pysam.AlignedSegment.fromstring(read, pysam.AlignmentHeader.from_dict(dict()))

        # Process and place our data on the output queue:
        out_queue.put(_segment_read(read, um))


def _reverse_segments(segments, read_length):
    """Reverses the segment order and coordinates to account for segmenting a reverse-complemented read."""
    rc_segments = []
    for seg in segments[::-1]:
        rc_segments.append(
            # Subtract 2 to account for inclusive coordinates on both ends:
            SegmentInfo(seg.name, read_length - seg.end - 2, read_length - seg.start - 2)
        )

    return rc_segments


def _collapse_annotations(path):
    """Collapses given path into a list of SegmentInfo objects."""
    last = ""
    start = 0
    segments = []
    for i, seg in enumerate(path):
        if seg != last:
            if i != 0:
                segments.append(SegmentInfo(last, start, i - 1))
            last = seg
            start = i
    # Don't forget the last one:
    segments.append(SegmentInfo(last, start, i - 1))
    return segments


def _segment_read(read, model):
    is_rc = False
    logp, ppath = annotate(model, read.query_sequence)

    rc_logp, rc_ppath = annotate(model, reverse_complement(read.query_sequence))
    if rc_logp > logp:
        logp = rc_logp
        ppath = rc_ppath
        is_rc = True
        logger.debug("Sequence scored better in RC: %s", read.query_name)

    return read.to_string(), ppath, logp, is_rc