import logging
import sys
import itertools
import re
import time
import collections
import math
import ctypes
import pickle
import numpy as np

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
__MP_BUFFER_LEN_BYTES = 5 * 1024 * 1024  # 5 Meg buffer
__MP_MESSAGE_DELIMITER = 0xDEADBEEF
__SENTINEL_VALUE = r"THE _MOST_ CROMULENT SENTINEL VALUE."

SEGMENTS_TAG = "SG"
SEGMENTS_RC_TAG = "RC"

logger = logging.getLogger(__name__)
click_log.basic_config(logger)

mod_name = "annotate"


# Create a new deque object that has put and get methods:
class Deque(collections.deque):

    _sleep_len = 0.000001

    def __init__(self, maxsize=None):
        super().__init__()
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


class ManyProducerSingleConsumerIpcQueue:
    """Shared memory backed queue class that can communicate between processes.
    Designed for many producers to be able to write to this queue, but only one consumer to read from it.
    Shared memory backing is implemented with multiprocessing.Array.
    """

    def __init__(self, size_bytes=1024*1024*5, sleep_time_s=0.0000001):
        self._buffer = mp.Array(size_or_initializer=size_bytes, typecode_or_type=ctypes.c_ubyte, lock=True)
        self._buffer_size = size_bytes

        # The following 3 counters are always updated together, so we use only 1 lock for them:
        self._num_messages_in_buffer = mp.Value(ctypes.c_ulong, 0, lock=True)
        self._cur_bytes_in_buffer = mp.Value(ctypes.c_ulong, 0, lock=False)
        self._unused_bytes_in_buffer = mp.Value(ctypes.c_ulong, size_bytes, lock=False)

        self._cur_buf_read_position = mp.Value(ctypes.c_ulong, 0, lock=True)
        self._cur_buf_write_position = mp.Value(ctypes.c_ulong, 0, lock=True)

        # we use a canonical delimiter here:
        # TODO: Update to use UUID lib for increased safety ala: uuid.uuid4().bytes
        self._delim = bytearray(b'\xDE\xAD\xBE\xEF')
        self._sleep_time_s = sleep_time_s

    def __len__(self):
        return self._num_messages_in_buffer

    def size(self):
        return self._buffer_size

    def put(self, obj):
        """Put an object onto this queue.
        The object must be able to be pickled.

        If the queue doesn't have enough free space to hold the object, blocks until enough space is available.
        If the queue is too small for the object, an exception is thrown."""

        data = pickle.dumps(obj, fix_imports=False)
        l_data = len(data) + len(self._delim)
        if l_data > self._buffer_size:
            raise BufferError(f"Given data cannot fit in buffer with delimiter: {l_data} > {self._buffer_size}")

        # Lock on the write position:
        with self._cur_buf_write_position.get_lock():
            # Here we peek at the bytes remaining in the buffer to wait for it to be small enough.
            # We can't write to this buffer while we're locked above, so we are guaranteed that the size will not grow.
            while l_data > self._unused_bytes_in_buffer.value:
                time.sleep(self._sleep_time_s)

            # Now we can write out our data to the buffer:
            with self._buffer.get_lock():
                # Do we have to write in multiple pieces because of the end of the buffer?
                if l_data > (self._buffer_size - self._cur_buf_write_position.value + 1):
                    # Write the first part of the data:
                    l_first_part = self._buffer_size - self._cur_buf_write_position.value + 1
                    self._buffer[self._cur_buf_write_position.value:
                                 self._cur_buf_write_position.value+l_first_part] = data[0:l_first_part]

                    # Write the second part:
                    l_second_part = len(data) - l_first_part
                    self._buffer[0:l_second_part] = data[l_first_part:]

                    # Write delimiter:
                    self._buffer[l_second_part:l_second_part+len(self._delim)] = self._delim[:]

                    # Update buffer info:
                    self._cur_buf_write_position.value = l_second_part + len(self._delim)
                else:
                    self._buffer[self._cur_buf_write_position.value:
                                 self._cur_buf_write_position.value+len(data)] = data[:]
                    self._cur_buf_write_position.value += len(data)
                    self._buffer[self._cur_buf_write_position.value:
                                 self._cur_buf_write_position.value+len(self._delim)] = self._delim[:]

                # Update our counts:
                with self._num_messages_in_buffer.get_lock():
                    self._unused_bytes_in_buffer.value -= l_data
                    self._cur_bytes_in_buffer.value += l_data
                    self._num_messages_in_buffer.value += 1

    def get(self):
        """Get an object off of this queue.
        If there are no objects on this queue, waits until one is available."""

        # Get our read lock:
        with self._cur_buf_read_position.get_lock():

            # Wait for a message to come in:
            while self._num_messages_in_buffer.value < 1:
                time.sleep(self._sleep_time_s)

            # Now that we have a message, we need to scan our buffer for the delimiter and grab the data:
            delim_start_pos = self._cur_buf_read_position.value
            delim_found = False
            while not delim_found:
                if delim_start_pos + len(self._delim) >= self._buffer_size:
                    # Delimiter could span the circular boundary.
                    # Get our data in two pieces then compare:
                    first_piece = self._buffer[delim_start_pos:]
                    second_piece = self._buffer[0:len(self._delim) - len(first_piece)]
                    if first_piece + second_piece == self._delim:
                        delim_found = True
                else:
                    if self._buffer[delim_start_pos:delim_start_pos+len(self._delim)] == self._delim:
                        # We found our delimiter.
                        delim_found = True
                delim_start_pos += len(self._delim)

            # Now delim_start_pos holds the delimiter start position.
            # It's possible that delim_start_pos is now before our current read position, which would mean we have
            # to read the data in 2 pieces:
            if delim_start_pos < self._cur_buf_read_position.value:
                # Get our data size:
                l_data = self._buffer_size - self._cur_buf_read_position.value + 1
                l_data += delim_start_pos

                data = bytearray(l_data)

                # Get the first part of the data:
                data[:self._buffer_size - self._cur_buf_read_position.value + 1] = \
                    self._buffer[self._cur_buf_read_position.value:]

                # Get the second part of the data:
                data[-delim_start_pos:] = self._buffer[:delim_start_pos]

            else:
                # Get the data up to the start of the delimiter:
                l_data = delim_start_pos - self._cur_buf_read_position.value
                data = bytearray(l_data)
                data[:l_data] = self._buffer[self._cur_buf_read_position.value:delim_start_pos]

            # Create our object:
            obj = pickle.loads(data, fix_imports=False)

            # Update our counts:
            self._cur_buf_read_position.value = delim_start_pos
            with self._num_messages_in_buffer.get_lock():
                self._unused_bytes_in_buffer.value += l_data
                self._cur_bytes_in_buffer.value -= l_data
                self._num_messages_in_buffer.value -= 1

        return obj


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

    for i in range(threads):
        if use_threads:
            p = threading.Thread(target=_worker_segmentation_fn, args=(input_data_queue, results, i))
        else:
            # Create an array to store our cross-process objects:
            if use_shared_memory:
                # TODO: Update to SharableList
                multi_process_arrays.append(
                    mp.Array(size_or_initializer=__MP_BUFFER_LEN_BYTES, typecode_or_type=ctypes.c_ubyte, lock=True)
                )
                ctypes.memset(multi_process_arrays[-1].get_obj(), 0x0, __MP_BUFFER_LEN_BYTES)
                p = mp.Process(target=_worker_segmentation_fn, args=(multi_process_arrays[i], results, i))

                # multi_process_arrays.append(ManyProducerSingleConsumerIpcQueue())
                # p = mp.Process(target=_worker_segmentation_fn, args=(multi_process_arrays[i], results, i))
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
        iter_data = itertools.chain(bam_file, (__SENTINEL_VALUE,) * threads)
        queue_put_time = 0
        cur_worker_thread = 0
        still_running = [True for _ in range(threads)]
        iter_st_time = time.time()
        for r in iter_data:
            # We have to adjust for our sentinel value if we've got to it:
            if r is not __SENTINEL_VALUE:
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
                                (r is __SENTINEL_VALUE and not still_running[worker_indx]):
                            continue
                        else:
                            raw = pickle.dumps(r, fix_imports=False)

                            with multi_process_arrays[worker_indx].get_lock():
                                multi_process_arrays[worker_indx][0:len(raw)] = raw[:]

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
        results.put(__SENTINEL_VALUE)
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


def __get_next_data_from_mp_circular_buffer(data_buffer,
                                            current_index,
                                            scratch_buffer=None,
                                            data_delim=__MP_MESSAGE_DELIMITER,
                                            read_chunk_size_bytes=1024):
    """Gets the next piece of data from the circular buffer."""
    found = False
    if not scratch_buffer:
        scratch_buffer = np.array(len(data_buffer), dtype=np.int8)
    scratch_indx = 0
    while not found:
        next_chunk_size = read_chunk_size_bytes
        if len(scratch_buffer) - scratch_indx - 1 < read_chunk_size_bytes:
            next_chunk_size = len(scratch_buffer) - scratch_indx - 1
        if len(data_buffer) - current_index - 1 < next_chunk_size:
            next_chunk_size = len(data_buffer) - current_index - 1

        with data_buffer.get_lock():
            scratch_buffer[scratch_indx:next_chunk_size] = data_buffer[current_index:next_chunk_size]


def __get_next_data_from_queue(q):
    try:
        # Multi-processing Array:
        first_byte = q[0]
        while first_byte == 0x0:
            time.sleep(__SLEEP_LEN_S)
            first_byte = q[0]
        with q.get_lock():
            raw_data = pickle.loads(q.get_obj(), fix_imports=False)
            ctypes.memset(q.get_obj(), 0x0, len(q))
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
            if raw_data == __SENTINEL_VALUE:
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
        if raw_data == __SENTINEL_VALUE:
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
