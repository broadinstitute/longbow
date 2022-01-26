import logging
import math
import sys
import itertools
import re
import time
import os
from collections import defaultdict

import click
import click_log
import tqdm

import ssw

import pysam
import multiprocessing as mp

import gzip
from construct import *

import longbow.utils.constants
from ..utils import bam_utils, barcode_utils
from ..utils.bam_utils import SegmentInfo
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("refine")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
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
    help="The model to use for annotation.  If not specified, it will be autodetected from "
         "the BAM header.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name "
         "and Longbow will attempt to read in the file and create a LibraryModel from it.  "
         "Longbow will assume the contents are the configuration of a LibraryModel as per "
         "LibraryModel.to_json()."
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
    "-b",
    "--barcode-tag",
    type=str,
    default=longbow.utils.constants.DEFAULT_MODEL,
    show_default=True,
    help="The barcode tag to adjust based on the allowlist"
)
@click.option(
    "-a",
    "--allow-list",
    type=click.Path(exists=True),
    help="list of allowed barcodes for specified tag (.txt, .txt.gz)",
)
@click.option(
    '-s',
    '--same-barcode-within-read',
    is_flag=True,
    default=False,
    show_default=True,
    help="Enforce constraint that all barcodes found within a single read are identical."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(threads, output_bam, model, force, barcode_tag, allow_list, same_barcode_within_read, input_bam):
    """Adjust model annotations based on a provided barcode allowlist."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Load barcode allow list
    barcodes = barcode_utils.load_barcode_allowlist(allow_list)

    # Configure process manager:
    # NOTE: We're using processes to overcome the Global Interpreter Lock.
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(threads)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_process_pool = []
    for _ in range(threads):
        p = mp.Process(
            target=_sub_process_work_fn, args=(process_input_data_queue, results, barcode_tag, barcodes, same_barcode_within_read)
        )
        p.start()
        worker_process_pool.append(p)

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        colour="green",
        file=sys.stderr,
        leave=False,
        disable=not sys.stdin.isatty(),
    ) as pbar:

        # Get our model:
        if model is None:
            lb_model = LibraryModel.from_json_obj(bam_utils.get_model_from_bam_header(bam_file.header))
        elif model is not None and LibraryModel.has_prebuilt_model(model):
            lb_model = LibraryModel.build_pre_configured_model(model)
        else:
            lb_model = LibraryModel.from_json_file(model)

        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[lb_model])

        # Start output worker:
        res = manager.dict({"num_reads_refined": 0, "num_reads": 0})
        output_worker = mp.Process(
            target=_refine_and_write_thread_fn,
            args=(
                results,
                out_header,
                output_bam,
                pbar,
                res,
                lb_model,
                barcode_tag,
                barcodes
            ),
        )
        output_worker.start()

        # Add in a sentinel value at the end of the queue - one for each subprocess - so we guarantee
        # that all subprocesses will exit:
        iter_data = itertools.chain(bam_file, (None,) * threads)
        for r in iter_data:
            if r is not None:
                process_input_data_queue.put(r.to_string())
            else:
                process_input_data_queue.put(r)

        # Wait for our input jobs to finish:
        for p in worker_process_pool:
            p.join()

        # Now that our input processes are done, we can add our exit sentinel onto the output queue and
        # wait for that process to end:
        results.put(None)
        output_worker.join()

    logger.info(
        f"Refined {res['num_reads_refined']} reads of {res['num_reads']} total."
    )
    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {res['num_reads']/(et - t_start):2.2f} reads/s.")


def get_segments(read):
    """Get the segments corresponding to a particular read by reading the segments tag information."""
    return read.to_string(), [
        SegmentInfo.from_tag(s) for s in read.get_tag(longbow.utils.constants.SEGMENTS_TAG).split(
            longbow.utils.constants.SEGMENT_TAG_DELIMITER)
    ]


def _refine_and_write_thread_fn(out_queue, out_bam_header, out_bam_file_name, pbar, res, model, barcode_tag, barcodes):
    """Thread / process fn to write out all our data."""

    with pysam.AlignmentFile(
        out_bam_file_name, "wb", header=out_bam_header
    ) as out_bam_file:
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
            read = raw_data

            # Condense the output annotations so we can write them out with indices:
            # segments = bam_utils.collapse_annotations(ppath)

            read = pysam.AlignedSegment.fromstring(read, out_bam_header)

            # Obligatory log message:
            # logger.debug(
            #     "Path for read %s (%2.2f)%s: %s",
            #     read.query_name,
            #     logp,
            #     " (RC)" if is_rc else "",
            #     segments,
            # )

            # Write our our read:
            # bam_utils.write_annotated_read(read, segments, is_rc, logp, model, ssw_aligner, out_bam_file)

            # Increment our counters:
            res["num_reads_refined"] += 1
            res["num_reads"] += 1

            pbar.update(1)


def _sub_process_work_fn(in_queue, out_queue, barcode_tag, barcodes, same_barcode_within_read):
    """Function to run in each subprocess.
    Extracts and returns all segments from an input read."""

    ssw_aligner = ssw.Aligner()

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
        read = pysam.AlignedSegment.fromstring(
            raw_data, pysam.AlignmentHeader.from_dict(dict())
        )

        if read.has_tag(longbow.utils.constants.SEGMENTS_TAG) and barcode_tag in read.get_tag(longbow.utils.constants.SEGMENTS_TAG):
            called_barcodes = defaultdict(int)

            segments = re.split(",", read.get_tag(longbow.utils.constants.SEGMENTS_TAG))
            for i, segment in enumerate(segments):
                if barcode_tag in segment:
                    tag, start, stop = re.split("[:-]", segment)
                    start = int(start)
                    stop = int(stop)

                    pad = 2 if read.get_tag("rq") > 0.0 else 4

                    bc = read.query_sequence[start - pad:stop + pad]

                    best_score = 0
                    best_alignment = None
                    best_barcode = None
                    for barcode in barcodes:
                        a = ssw_aligner.align(query=bc, reference=barcode, revcomp=False)

                        if a.score > best_score:
                            best_score = a.score
                            best_alignment = a
                            best_barcode = barcode

                    print(f"{i} {start} {stop} {barcode_tag} {segment} {read.get_tag('rq')} {best_barcode}")
                    print(best_alignment.alignment_report())
                    print("")

                    called_barcodes[best_barcode] += 1

            print(f"{read.query_name} {called_barcodes}")
            logger.info("")

        # Process and place our data on the output queue:
        out_queue.put(read.to_string())


def _refine_read(read, model):
    is_rc = False
    logp, ppath = model.annotate(read.query_sequence)

    rc_logp, rc_ppath = model.annotate(bam_utils.reverse_complement(read.query_sequence))

    # print(f"Forward Path: {logp}: {ppath}")
    # print(f"Reverse Path: {rc_logp}: {rc_ppath}")

    if rc_logp > logp:
        logp = rc_logp
        ppath = rc_ppath
        is_rc = True
        logger.debug("Sequence scored better in RC: %s", read.query_name)

    return read.to_string(), ppath, logp, is_rc
