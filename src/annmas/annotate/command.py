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


class ManyProducerSingleConsumerIpcQueue:
    """Shared memory backed queue class that can communicate between processes.
    Designed for many producers to be able to write to this queue, but only one consumer to read from it.
    Shared memory backing is implemented with multiprocessing.Array.
    """

    __instance_counter = 0

    def __init__(self, size_bytes=1024*1024*5, sleep_time_s=0.0000001):
        self._buffer = mp.Array(size_or_initializer=size_bytes, typecode_or_type=ctypes.c_ubyte, lock=True)
        ctypes.memset(self._buffer.get_obj(), 0x0, size_bytes)
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

        self._label = ManyProducerSingleConsumerIpcQueue.__instance_counter
        ManyProducerSingleConsumerIpcQueue.__instance_counter += 1

        # So we can track how long it takes for us to pickle/unpickle our data:
        self._total_pickle_time_s = mp.Value(ctypes.c_double, 0, lock=True)
        self._total_unpickle_time_s = mp.Value(ctypes.c_double, 0, lock=True)

        # To track sleep time:
        self._total_time_spent_waiting_in_put_s = mp.Value(ctypes.c_double, 0, lock=True)
        self._total_time_spent_waiting_in_get_s = mp.Value(ctypes.c_double, 0, lock=True)

    def __len__(self):
        return self._num_messages_in_buffer.value

    def total_pickle_time(self):
        return self._total_pickle_time_s.value

    def total_unpickle_time(self):
        return self._total_unpickle_time_s.value

    def total_put_wait_time(self):
        return self._total_time_spent_waiting_in_put_s.value

    def total_get_wait_time(self):
        return self._total_time_spent_waiting_in_get_s.value

    def size(self):
        return self._buffer_size

    def put(self, obj):
        """Put an object onto this queue.
        The object must be able to be pickled.

        If the queue doesn't have enough free space to hold the object, blocks until enough space is available.
        If the queue is too small for the object, an exception is thrown."""

        st = time.time()
        data = pickle.dumps(obj, fix_imports=False)
        et = time.time()
        self._total_pickle_time_s.acquire()
        self._total_pickle_time_s.value += et - st
        self._total_pickle_time_s.release()

        packet = data + self._delim
        l_packet = len(data) + len(self._delim)
        if l_packet > self._buffer_size:
            raise BufferError(f"Given data cannot fit in buffer with delimiter: {l_packet} > {self._buffer_size}")

        # Lock on the write position:
        try:
            logger.debug(f"IPCQ {self._label} about to acquire write position lock.")
            self._cur_buf_write_position.acquire()
            logger.debug(f"IPCQ {self._label} acquired write position lock.")

            # Here we peek at the bytes remaining in the buffer to wait for it to be small enough.
            # We can't write to this buffer while we're locked above, so we are guaranteed that the size will not grow.
            have_logged = False
            st = time.time()
            while l_packet > self._unused_bytes_in_buffer.value:
                if not have_logged:
                    logger.debug("Waiting for space to place object (%d > %d/%d free bytes, items in queue: %d).",
                                 l_packet, self._unused_bytes_in_buffer.value, self._buffer_size,
                                 self._num_messages_in_buffer.value)
                    have_logged = True
                time.sleep(self._sleep_time_s)
            et = time.time()
            self._total_time_spent_waiting_in_put_s.acquire()
            self._total_time_spent_waiting_in_put_s.value += et - st
            self._total_time_spent_waiting_in_put_s.release()

            # Now we can write out our data to the buffer:
            # with self._buffer.get_lock():

            # Do we have to write in multiple pieces because of the end of the buffer?
            if l_packet > (self._buffer_size - self._cur_buf_write_position.value):
                # Write the first part of the data:
                l_first_part = self._buffer_size - self._cur_buf_write_position.value

                self._buffer[self._cur_buf_write_position.value:
                             self._cur_buf_write_position.value+l_first_part] = packet[:l_first_part]

                # Write the second part:
                l_second_part = l_packet - l_first_part
                self._buffer[:l_second_part] = packet[l_first_part:]

                # Update buffer info:
                self._cur_buf_write_position.value = l_second_part
            else:
                self._buffer[self._cur_buf_write_position.value:
                             self._cur_buf_write_position.value+l_packet] = packet
                self._cur_buf_write_position.value += l_packet

            # Update our counts:
            try:
                logger.debug(f"IPCQ {self._label} about to acquire num message lock (put).")
                self._num_messages_in_buffer.acquire()
                logger.debug(f"IPCQ {self._label} acquired num message lock (put).")

                self._unused_bytes_in_buffer.value -= l_packet
                self._cur_bytes_in_buffer.value += l_packet
                self._num_messages_in_buffer.value += 1
            finally:
                logger.debug(f"IPCQ {self._label} about to release num message lock (put).")
                self._num_messages_in_buffer.release()
                logger.debug(f"IPCQ {self._label} released num message lock (put).")
        finally:
            self._cur_buf_write_position.release()

    def get(self):
        """Get an object off of this queue.
        If there are no objects on this queue, waits until one is available."""

        # Wait for a message to come in:
        st = time.time()
        while self._num_messages_in_buffer.value < 1:
            time.sleep(self._sleep_time_s)
        et = time.time()
        self._total_time_spent_waiting_in_get_s.acquire()
        self._total_time_spent_waiting_in_get_s.value += et - st
        self._total_time_spent_waiting_in_get_s.release()

        # Get our read lock:
        try:
            logger.debug(f"IPCQ {self._label} about to acquire read position lock.")
            self._cur_buf_read_position.acquire()
            logger.debug(f"IPCQ {self._label} acquired read position lock.")

            # Now that we have a message, we need to scan our buffer for the delimiter and grab the data:
            delim_start_pos = self._cur_buf_read_position.value
            delim_found = False
            while not delim_found:
                if delim_start_pos + len(self._delim) > self._buffer_size:
                    # Delimiter could span the circular boundary.
                    # Get our data in two pieces then compare:
                    first_piece = self._buffer[delim_start_pos:]
                    second_piece = self._buffer[0:len(self._delim) - len(first_piece)]
                    if self._are_arrays_equal(first_piece + second_piece, self._delim):
                        delim_found = True
                else:
                    if self._are_arrays_equal(self._buffer[delim_start_pos:delim_start_pos + len(self._delim)],
                                              self._delim):
                        # We found our delimiter.
                        delim_found = True
                delim_start_pos = delim_start_pos + 1 if delim_start_pos < self._buffer_size else 0

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
            st = time.time()
            obj = pickle.loads(data, fix_imports=False)
            et = time.time()
            self._total_unpickle_time_s.acquire()
            self._total_unpickle_time_s.value += et - st
            self._total_unpickle_time_s.release()

            # Update our counts:
            self._cur_buf_read_position.value = delim_start_pos + len(self._delim) - 1
            try:
                logger.debug(f"IPCQ {self._label} about to acquire num message lock (get).")
                self._num_messages_in_buffer.acquire()
                logger.debug(f"IPCQ {self._label} acquired num message lock (get).")

                self._unused_bytes_in_buffer.value += l_data + len(self._delim) - 1
                self._cur_bytes_in_buffer.value -= l_data + len(self._delim) - 1
                self._num_messages_in_buffer.value -= 1
            finally:
                logger.debug(f"IPCQ {self._label} about to release num message lock (get).")
                self._num_messages_in_buffer.release()
                logger.debug(f"IPCQ {self._label} released num message lock (get).")

        finally:
            logger.debug(f"IPCQ {self._label} about to release read position lock.")
            self._cur_buf_read_position.release()
            logger.debug(f"IPCQ {self._label} released read position lock.")

        return obj

    @staticmethod
    def _are_arrays_equal(it1, it2):
        if len(it1) != len(it2):
            return False
        for i in range(len(it1)):
            if it1[i] != it2[i]:
                return False
        return True


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

    t_start = time.time()
    logger.info(f"annmas: {mod_name} started")

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    m = build_default_model()
    if model is not None:
        m.from_yaml(model)
        logger.info(f"Using pretrained annotation model {model}")
    else:
        logger.info("Using default annotation model")

    # Create queues for data:
    queue_size = threads*2 if threads < 10 else 20
    manager = mp.Manager()
    input_data_queue = manager.Queue(maxsize=queue_size)
    results = manager.Queue()

    # input_data_queue = ManyProducerSingleConsumerIpcQueue()
    # results = ManyProducerSingleConsumerIpcQueue()

    # Start worker sub-processes:
    worker_pool = []
    # input_data_queue = ManyProducerSingleConsumerIpcQueue()

    for i in range(threads):
        p = mp.Process(target=_worker_segmentation_fn, args=(input_data_queue, results, i, m))
        p.start()
        worker_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(desc="Progress", unit=" read", colour="green", file=sys.stdout) as pbar:

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
        output_worker = mp.Process(target=_write_thread_fn, args=(results, out_header, output_bam, pbar))
        output_worker.start()

        # Add in a `None` sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        iter_data = itertools.chain(bam_file, (__SENTINEL_VALUE,) * threads)
        for r in iter_data:
            # We have to adjust for our sentinel value if we've got to it:
            if r is not __SENTINEL_VALUE:
                r = r.to_string()
            input_data_queue.put(r)

        logger.info("Finished reading data and sending it to sub-processes.")
        logger.info("Waiting for sub-processes to finish...")

        # Wait for our input jobs to finish:
        for p in worker_pool:
            p.join()

        logger.debug("All workers stopped.")
        logger.debug("Terminating output process.")

        # Now that our input processes are done, we can add our exit sentinel onto the output queue and
        # wait for that process to end:
        results.put(__SENTINEL_VALUE)
        output_worker.join()

    logger.info(f"annmas: {mod_name} finished.")
    logger.info(f"annmas: elapsed time: %2.5fs.", time.time() - t_start)


def _write_thread_fn(out_queue, out_bam_header, out_bam_file_name, pbar):
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
            if raw_data == __SENTINEL_VALUE:
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
        if raw_data == __SENTINEL_VALUE:
            break
        # Should really never be None, but just in case:
        elif raw_data is None:
            continue

        # Unpack our data here:
        read = raw_data
        read = pysam.AlignedSegment.fromstring(read, pysam.AlignmentHeader.from_dict(dict()))

        # Process and place our data on the output queue:
        segment_info = _segment_read(read, model)

        out_queue.put(segment_info)
        num_reads_segmented += 1

    logger.info(f"Worker %d: Num reads segmented: %d", worker_num, num_reads_segmented)


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
