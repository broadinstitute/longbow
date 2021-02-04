import logging
import sys
import itertools
import re
import time
import collections
import math
import ctypes
import pickle

import queue
from inspect import getframeinfo, currentframe, getdoc

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp
import threading

from ..utils.model import build_default_model
from ..utils.model import reverse_complement
from ..utils.model import annotate

import annmas.fmq as fmq

from ..meta import VERSION

__SLEEP_LEN_S = 0.000001
SEGMENTS_TAG = "SG"
SEGMENTS_RC_TAG = "RC"

_SENTINEL_VALUE = r"THE _MOST_ CROMULENT SENTINEL VALUE."

logger = logging.getLogger(__name__)
click_log.basic_config(logger)

mod_name = "annotate"


# Create a new deque object that has put and get methods:
class Deque(collections.deque):

    _sleep_len = 0.000001

    def __init__(self, maxsize=None):
        self.maxsize = maxsize if maxsize else math.inf

    def put(self, o):
        while len(self) >= self.maxsize:
            # Sleep for arbitrarily small amount of time:
            time.sleep(Deque._sleep_len)
        self.append(o)

    def get(self):
        while len(self) == 0:
            # Sleep for arbitrarily small amount of time to avoid stack trace:
            time.sleep(Deque._sleep_len)
        return self.popleft()


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
    """Annotate reads in a BAM file with segments from the model."""

    use_threads = False
    use_shared_memory = True

    t_start = time.time()
    logger.info(f"annmas: {mod_name} started")

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    if use_threads:
        logger.info(f"Running with {threads} subthread(s)")
    else:
        logger.info(f"Running with {threads} subprocess(es)")

    m = build_default_model()
    if model is not None:
        m.from_yaml(model)
        logger.info(f"Using pretrained annotation model {model}")
    else:
        logger.info("Using default annotation model")

    # Create queues for data:
    queue_size = threads*2 if threads < 10 else 20
    if use_threads:
        input_data_queue = queue.Queue(maxsize=queue_size)
        results = queue.Queue()
        # input_data_queue = Deque()
        # results = Deque()
    else:
        manager = mp.Manager()
        input_data_queue = manager.Queue(maxsize=queue_size)
        results = manager.Queue()

    # Start worker sub-processes:
    worker_pool = []
    multi_process_arrays = []
    mp_array_len = 1024 * 1024
    for i in range(threads):
        if use_threads:
            p = threading.Thread(target=_worker_segmentation_fn, args=(input_data_queue, results, i))
        else:
            # Create an array to store our cross-process objects:
            if use_shared_memory:
                # TODO: Update to SharableList
                multi_process_arrays.append(
                    mp.Array(size_or_initializer=mp_array_len, typecode_or_type=ctypes.c_ubyte, lock=True)
                )
                ctypes.memset(multi_process_arrays[-1].get_obj(), 0x0, mp_array_len)
                p = mp.Process(target=_worker_segmentation_fn, args=(multi_process_arrays[i], results, i))
            else:
                p = mp.Process(target=_worker_segmentation_fn, args=(input_data_queue, results, i))
        p.start()
        worker_pool.append(p)

    serialize_time_s = 0

    pbar = None

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file:  #, tqdm.tqdm(desc="Progress", unit=" read", colour="green", file=sys.stdout) as pbar:

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
        if use_threads:
            output_worker = threading.Thread(target=_write_thread_fn, args=(results, out_header, output_bam, pbar))
            output_worker.start()
        else:
            output_worker = mp.Process(target=_write_thread_fn, args=(results, out_header, output_bam, pbar))
            output_worker.start()

        # Add in a `None` sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        iter_data = itertools.chain(bam_file, (_SENTINEL_VALUE,) * threads)
        queue_put_time = 0
        cur_worker_thread = 0
        still_running = [True for _ in range(threads)]
        iter_st_time = time.time()
        for r in iter_data:
            # We have to adjust for our sentinel value if we've got to it:
            if r is not _SENTINEL_VALUE:
                st = time.time()
                rs = r.to_string()
                serialize_time_s += time.time() - st
                r = (rs, m)

            st = time.time()
            if use_shared_memory and not use_threads:
                placed = False
                while not placed:
                    for i in range(threads):
                        worker_indx = (cur_worker_thread + i) % threads
                        if multi_process_arrays[worker_indx][0] != 0x0 or \
                                (r is _SENTINEL_VALUE and not still_running[worker_indx]):
                            continue
                        else:
                            raw = pickle.dumps(r, fix_imports=False)

                            multi_process_arrays[worker_indx].acquire()
                            multi_process_arrays[worker_indx][0:len(raw)] = raw[:]
                            multi_process_arrays[worker_indx].release()

                            cur_worker_thread = (worker_indx + 1) % threads
                            placed = True
                            # if r is _SENTINEL_VALUE:
                            #     still_running[worker_indx] = False
                            break
                    if not placed:
                        time.sleep(__SLEEP_LEN_S)
            else:
                input_data_queue.put((rs, m))
            queue_put_time += time.time() - st
            # print("Still running: " + str(still_running))

        iter_time = time.time() - iter_st_time
        logger.info(f"Queue put time: %2.5f", queue_put_time)
        logger.info(f"Iteration time: %2.5f", iter_time)

        # Wait for our input jobs to finish:
        for p in worker_pool:
            p.join()

        logger.debug("All workers stopped.")

        # Now that our input processes are done, we can add our exit sentinel onto the output queue and
        # wait for that process to end:
        results.put(_SENTINEL_VALUE)
        output_worker.join()

    logger.info(f"annmas: {mod_name} finished.")
    logger.info(f"annmas: elapsed time: %2.5fs.", time.time() - t_start)
    logger.info("")
    logger.info(f"Read serialize to string time: %2.5f", serialize_time_s)
# def main(model, threads, output_bam, input_bam):
#     """Annotate reads in a BAM file with segments from the model."""
#
#     # This is the single-threaded version:
#
#     t_start = time.time()
#     logger.info(f"annmas: {mod_name} started")
#
#     m = build_default_model()
#     if model is not None:
#         m.from_yaml(model)
#         logger.info(f"Using pretrained annotation model {model}")
#     else:
#         logger.info("Using default annotation model")
#
#     segmentation_time_s = 0
#
#     pysam.set_verbosity(0)  # silence message about the .bai file not being found
#     with pysam.AlignmentFile(
#         input_bam, "rb", check_sq=False, require_index=False
#     ) as bam_file, tqdm.tqdm(desc="Progress", unit=" read", colour="green", file=sys.stdout) as pbar:
#
#         # Get our header from the input bam file:
#         out_bam_header_dict = bam_file.header.to_dict()
#
#         # Add our program group to it:
#         pg_dict = {
#             "ID": f"annmas-annotate-{VERSION}",
#             "PN": "annmas",
#             "VN": f"{VERSION}",
#             # Use reflection to get the doc string for this main function for our header:
#             "DS": getdoc(globals()[getframeinfo(currentframe()).function]),
#             "CL": " ".join(sys.argv),
#         }
#         if "PG" in out_bam_header_dict:
#             out_bam_header_dict["PG"].append(pg_dict)
#         else:
#             out_bam_header_dict["PG"] = [pg_dict]
#         out_header = pysam.AlignmentHeader.from_dict(out_bam_header_dict)
#
#         with pysam.AlignmentFile(
#                 output_bam, "wb", header=out_header
#         ) as out_bam_file:
#             iter_st_time = time.time()
#             for read in bam_file:
#                 st = time.time()
#                 segment_info = _segment_read(read, m)
#                 segmentation_time_s += time.time() - st
#
#                 read, ppath, logp, is_rc = segment_info
#                 read = pysam.AlignedSegment.fromstring(read, out_header)
#
#                 # Condense the output annotations so we can write them out with indices:
#                 segments = _collapse_annotations(ppath)
#
#                 # Obligatory log message:
#                 logger.debug(
#                     "Path for read %s (%2.2f)%s: %s",
#                     read.query_name,
#                     logp,
#                     " (RC)" if is_rc else "",
#                     segments,
#                 )
#
#                 # Set our tag and write out the read to the annotated file:
#                 read.set_tag(SEGMENTS_TAG, "|".join([s.to_tag() for s in segments]))
#
#                 # If we're reverse complemented, we make it easy and just reverse complement the read and add a tag
#                 # saying that the read was RC:
#                 read.set_tag(SEGMENTS_RC_TAG, is_rc)
#                 if is_rc:
#                     quals = read.query_qualities[::-1]
#                     seq = reverse_complement(read.query_sequence)
#                     read.query_sequence = seq
#                     read.query_qualities = quals
#                 out_bam_file.write(read)
#
#         iter_time = time.time() - iter_st_time
#         logger.info(f"Iteration time: %2.5f", iter_time)
#
#     logger.info(f"annmas: {mod_name} finished.")
#     logger.info(f"annmas: segmentation time: %2.5fs.", segmentation_time_s)
#     logger.info(f"annmas: elapsed time: %2.5fs.", time.time() - t_start)


def __get_next_data_from_queue(q):
    try:
        # Multi-processing Array:
        first_byte = q[0]
        while first_byte == 0x0:
            time.sleep(__SLEEP_LEN_S)
            first_byte = q[0]
        q.acquire()
        raw_data = pickle.loads(q.get_obj(), fix_imports=False)
        ctypes.memset(q.get_obj(), 0x0, len(q))
        q.release()
    # this happens if you can't subscript the data, which you can't with queues:
    except TypeError:
        raw_data = q.get()

    return raw_data


def _write_thread_fn(out_queue, out_bam_header, out_bam_file_name, pbar):
    """Thread / process fn to write out all our data."""

    num_reads_annotated = 0
    num_sections = 0

    deserialize_time_s = 0
    wait_time_s = 0

    with pysam.AlignmentFile(
            out_bam_file_name, "wb", header=out_bam_header
    ) as out_bam_file:

        while True:
            # Wait for some output data:
            st = time.time()
            raw_data = __get_next_data_from_queue(out_queue)
            wait_time_s += time.time() - st

            # Check for exit sentinel:
            if raw_data == _SENTINEL_VALUE:
                break
            # Should really never be None, but just in case:
            elif raw_data is None:
                continue

            # Unpack data:
            read, ppath, logp, is_rc = raw_data
            st = time.time()
            read = pysam.AlignedSegment.fromstring(read, out_bam_header)
            deserialize_time_s += time.time() - st

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

            # pbar.update(1)

    logger.info(
        f"annmas {mod_name}: annotated {num_reads_annotated} reads with {num_sections} total sections."
    )
    logger.info(f"Read from string deserialize time: %2.5f", deserialize_time_s)
    logger.info(f"Output process wait on input time: %2.5f", wait_time_s)
    return


def _worker_segmentation_fn(in_queue, out_queue, worker_num):
    """Function to run in each subthread / subprocess.
    Segments each read and place the segments in the output queue."""

    wait_time_s = 0
    segmentation_time_s = 0
    num_reads_segmented = 0

    while True:
        # Wait until we get some data.
        # Note: Because we have a sentinel value None inserted at the end of the input data for each
        #       subprocess, we don't have to add a timeout - we're guaranteed each process will always have
        #       at least one element.
        st = time.time()
        raw_data = __get_next_data_from_queue(in_queue)
        wait_time_s += time.time() - st

        # Check for exit sentinel:
        if raw_data == _SENTINEL_VALUE:
            break
        # Should really never be None, but just in case:
        elif raw_data is None:
            continue

        # Unpack our data here:
        try:
            read, um = raw_data
        except TypeError as e:
            print(raw_data)
            raise e
        read = pysam.AlignedSegment.fromstring(read, pysam.AlignmentHeader.from_dict(dict()))

        # Process and place our data on the output queue:
        st = time.time()
        segment_info = _segment_read(read, um)
        segmentation_time_s += time.time() - st

        out_queue.put(segment_info)
        num_reads_segmented += 1

    logger.info(f"Worker %d: Num reads segmented: %d", worker_num, num_reads_segmented)
    logger.info(f"Worker %d: Total segmentation time: %2.5f", worker_num, segmentation_time_s)
    logger.info(f"Worker %d: Input process wait on input time: %2.5f", worker_num, wait_time_s)


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
