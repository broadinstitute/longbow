import logging
import sys
import itertools
import re
import time
import collections
import os

from inspect import getframeinfo, currentframe, getdoc

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp

import gzip
from construct import *

from ..utils.model import build_default_model
from ..utils.model import reverse_complement
from ..utils.model import annotate

from ..meta import VERSION

SEGMENTS_TAG = "SG"
SEGMENTS_RC_TAG = "RC"

__SEGMENT_TAG_DELIMITER = ","

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("annotate")
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
    "-m",
    "--model",
    required=False,
    type=click.Path(exists=True),
    help="pre-trained model to apply",
)
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
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(model, pbi, threads, output_bam, input_bam):
    """Annotate reads in a BAM file with segments from the model."""

    t_start = time.time()

    logger.info("Invoked via: annmas %s", " ".join(sys.argv[1:]))

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    m = build_default_model()
    if model is not None:
        m.from_yaml(model)
        logger.info(f"Using pretrained annotation model {model}")
    else:
        logger.info("Using default annotation model")

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = _load_read_count(pbi)
        logger.info("Annotating %d reads", read_count)

    # Create queues for data:
    queue_size = threads * 2 if threads < 10 else 20
    manager = mp.Manager()
    input_data_queue = manager.Queue(maxsize=queue_size)
    results = manager.Queue()

    # input_data_queue = ManyProducerSingleConsumerIpcQueue()
    # results = ManyProducerSingleConsumerIpcQueue()

    # Start worker sub-processes:
    worker_pool = []
    # input_data_queue = ManyProducerSingleConsumerIpcQueue()

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
        out_bam_header_dict = bam_file.header.to_dict()

        # Add our program group to it:
        pg_dict = {
            "ID": f"annmas-annotate-{VERSION}",
            "PN": "annmas",
            "VN": f"{VERSION}",
            # Use reflection to get the doc string for this main function for our header:
            "DS": getdoc(globals()[getframeinfo(currentframe()).function]),
            "CL": " ".join(sys.argv),
        }
        if "PG" in out_bam_header_dict:
            out_bam_header_dict["PG"].append(pg_dict)
        else:
            out_bam_header_dict["PG"] = [pg_dict]
        out_header = pysam.AlignmentHeader.from_dict(out_bam_header_dict)

        # Start output worker:
        res = manager.dict({"num_reads_annotated": 0, "num_sections": 0})
        output_worker = mp.Process(
            target=_write_thread_fn, args=(results, out_header, output_bam, not sys.stdin.isatty(), res, read_count)
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


def _write_thread_fn(out_queue, out_bam_header, out_bam_file_name, disable_pbar, res, read_count):
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
            read.set_tag(SEGMENTS_TAG, __SEGMENT_TAG_DELIMITER.join([s.to_tag() for s in segments]))

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


def _reverse_segments(segments, read_length):
    """Reverses the segment order and coordinates to account for segmenting a reverse-complemented read."""
    rc_segments = []
    for seg in segments[::-1]:
        rc_segments.append(
            # Subtract 2 to account for inclusive coordinates on both ends:
            SegmentInfo(
                seg.name, read_length - seg.end - 2, read_length - seg.start - 2
            )
        )

    return rc_segments


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


def _load_read_count(pbi_file):
    """
    Compute file offsets for specified read names
    """

    # Decode PacBio .pbi file.  This is not a full decode of the index, only the parts we need
    # until we get to the read count.
    # More on index format at https://pacbiofileformats.readthedocs.io/en/9.0/PacBioBamIndex.html .

    fmt = Struct(
        # Header
        "magic" / Const(b"PBI\x01"),
        "version_patch" / Int8ul,
        "version_minor" / Int8ul,
        "version_major" / Int8ul,
        "version_empty" / Int8ul,
        "pbi_flags" / Int16ul,
        "n_reads" / Int32ul,
        )

    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        return idx_contents.n_reads

def _get_segments(read):
    """Get the segments corresponding to a particular read by reading the segments tag information."""
    return read.to_string(), [
        SegmentInfo.from_tag(s) for s in read.get_tag(SEGMENTS_TAG).split(__SEGMENT_TAG_DELIMITER)
    ]
