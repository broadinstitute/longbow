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
from ..utils.cli_utils import format_obnoxious_warning_message


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("pad")
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
    help="The barcode tag to adjust"
)
@click.option(
    "-n",
    "--new-barcode-tag",
    type=str,
    default=longbow.utils.constants.READ_RAW_UMI_TAG,
    show_default=True,
    help="The barcode tag into which to put the adjusted value."
)
@click.option(
    "-e",
    "--expand",
    type=int,
    default=2,
    help="Expand tag by specified number of bases",
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(threads, output_bam, model, force, barcode_tag, new_barcode_tag, expand, input_bam):
    """Pad tag by specified number of adjacent bases from the read."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    logger.info(f"Expanding tag {barcode_tag} by {expand} bases")
    logger.info(f"Writing expanded tag to: {new_barcode_tag}")

    # Configure process manager:
    # NOTE: We're using processes to overcome the Global Interpreter Lock.
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(threads)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_process_pool = []
    for _ in range(threads):
        p = mp.Process(
            target=_pass_read_fn, args=(process_input_data_queue, results)
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

        # Verify that the given model actually has the barcode to change:
        if not lb_model.has_annotation_tag(barcode_tag):
            print(f"ERROR: Could not determine {lb_model.name} model segment from tag name: {barcode_tag}",
                  file=sys.stderr)
            sys.exit(1)

        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[lb_model])

        # Start output worker:
        res = manager.dict({"num_reads_refined": 0, "num_reads": 0})
        output_worker = mp.Process(
            target=_expand_tag_fn,
            args=(
                results,
                out_header,
                output_bam,
                pbar,
                res,
                lb_model,
                barcode_tag,
                new_barcode_tag,
                expand
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

    logger.info( f"Refined {res['num_reads_refined']} reads of {res['num_reads']} total.")
    
    et = time.time()
    logger.info(f"Done. Elapsed time: {et - t_start:2.2f}s. "
                f"Overall processing rate: {res['num_reads']/(et - t_start):2.2f} reads/s.")

    if res['num_reads_refined'] == 0:
        logger.warning(format_obnoxious_warning_message("No reads were refined / padded.  This is very likely a misconfiguration."))


def _expand_tag_fn(out_queue, out_bam_header, out_bam_file_name, pbar, res, lb_model,
                   barcode_tag, new_barcode_tag, expand):
    """Thread / process fn to expand a tag and write out all our data."""

    # Get the segment to pad out:
    # NOTE: by now we know this tag is in the model:
    barcode_seg_name = lb_model.get_segment_name_for_annotation_tag(barcode_tag)

    with pysam.AlignmentFile(
        out_bam_file_name, "wb", header=out_bam_header
    ) as out_bam_file:
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
            read = pysam.AlignedSegment.fromstring(read, out_bam_header)

            read_is_refined = False

            if read.has_tag(barcode_tag):
                segments = read.get_tag(longbow.utils.constants.SEGMENTS_TAG).split(",")

                for segment in segments:
                    if segment.startswith(barcode_seg_name + ":"):
                        tag, start, stop = re.split("[:-]", segment)

                        tag_bc = read.get_tag(barcode_tag)

                        for query_sequence in [read.query_sequence, bam_utils.reverse_complement(read.query_sequence)]:
                            # Must add 1 to include the final base:
                            pos = query_sequence.find(read.get_tag(barcode_tag), int(start), int(stop)+1)

                            if pos >= 0:
                                start_pos = pos - expand
                                if start_pos < 0:
                                    start_pos = 0
                                end_pos = pos+len(tag_bc)+expand

                                new_bc = query_sequence[start_pos:end_pos]

                                # old_bc = query_sequence[pos:pos + len(tag_bc)]
                                # logger.info(f"{read.query_name} {pos}")
                                # logger.info(f" - {expand * ' '}{tag_bc}")
                                # logger.info(f" - {expand * ' '}{old_bc}")
                                # logger.info(f" - {new_bc}")
                                # logger.info("")

                                read.set_tag(new_barcode_tag, new_bc)
                                read_is_refined = True

                                # We only need to look in the correct direction of the read.
                                # No need to continue here:
                                break

            # Write our our read:
            out_bam_file.write(read)

            # Increment our counters:
            res["num_reads"] += 1
            if read_is_refined > 0:
                res["num_reads_refined"] += 1 if read_is_refined else 0

            pbar.update(1)


def _pass_read_fn(in_queue, out_queue):
    """Function to run in each subprocess.
    Just pass reads to the padding thread."""

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
        # read = pysam.AlignedSegment.fromstring(
        #     raw_data, pysam.AlignmentHeader.from_dict(dict())
        # )

        # if read.has_tag(longbow.utils.constants.SEGMENTS_TAG) and barcode_tag in read.get_tag(longbow.utils.constants.SEGMENTS_TAG):
        #     segments = read.get_tag(longbow.utils.constants.SEGMENTS_TAG).split(",")
        #     num_segments = len(segments)

        #     for i, segment in enumerate(segments):
        #         if barcode_tag in segment:
        #             tag, start, stop = re.split("[:-]", segment)
        #             start = int(start)
        #             stop = int(stop)

        #             old_bc = read.query_sequence[start:stop]
        #             new_bc = read.query_sequence[start - expand:stop + expand]

        #             logger.info(f'{tag} {old_bc} {new_bc}')
        #             logger.info("")

        # # Process and place our data on the output queue:
        # out_queue.put((read.to_string(), refined_segments, num_segments))

        # out_queue.put(read.to_string())
        out_queue.put(raw_data)
