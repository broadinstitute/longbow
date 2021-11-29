import logging
import sys
import itertools

import time

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp

from inspect import getframeinfo, currentframe, getdoc

import longbow.utils.constants
from ..utils import bam_utils
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel

from ..annotate.command import get_segments

from ..meta import VERSION

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("scsplit")
click_log.basic_config(logger)


__DEFAULT_DUMMY_CBC = "CTGCCTAACCTGATCC"
__DEFAULT_OUT_BASE_NAME = logger.name
__DEFAULT_UMI_LENGTH = 10


__OUT_READ_FILE_SUFFIX = "_mates"
__OUT_WHITELIST_FILE_SUFFIX = "_whitelist.txt"


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
    "--output-base-name",
    default=__DEFAULT_OUT_BASE_NAME,
    type=str,
    help=f"base name for output files [default: {__DEFAULT_OUT_BASE_NAME}]",
)
@click.option(
    "-c",
    "--cell-barcode",
    default=__DEFAULT_DUMMY_CBC,
    type=str,
    help=f"dummy cell barcode to use for the dataset [default: {__DEFAULT_DUMMY_CBC}, "
         f"length: {len(__DEFAULT_DUMMY_CBC)}]",
)
@click.option(
    "-u",
    "--umi-length",
    default=__DEFAULT_UMI_LENGTH,
    type=int,
    show_default=True,
    help=f"length of the UMI from this library",
)
@click.option(
    "-b",
    "--write-bam",
    is_flag=True,
    default=False,
    show_default=True,
    help=f"Write out an annotated bam file in addition to the mates files.",
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force scsplit to run on the input bam without checking for compatibility."
)
@click.option(
    "-m",
    "--model",
    default=longbow.utils.constants.DEFAULT_MODEL,
    show_default=True,
    help="The model to use for annotation.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(threads, output_base_name, cell_barcode, umi_length, force, model, write_bam, input_bam):
    """Create files for use in `alevin` for single-cell analysis.
    This tool coerces a set of reads from a single source into a format that `alevin` can ingest.
    
    Segment names are assumed to be those in the default model (utils/model.py).
    
    INPUT_BAM should contain reads that have been processed by `longbow segment`.

    The output from this tool consists of several files:
        OUTPUT_BASE_NAME_mates_1.fastq:
            A file containing partial sequences for all reads in the given input file.  These partial reads consist of the 
            dummy cell barcode + detected UMI for each read in the given input file.
        OUTPUT_BASE_NAME_mates_2.fastq:
            A file containing partial sequences for all reads in the given input file.  These partial reads consist of the 
            transcript sequences for all reads in the given input file.  Transcript sequences include data after the UMI 
            and before the Poly-A tail.  All bases outside of this range are excluded from the output.
        OUTPUT_BASE_NAME_whitelist.txt:
            A whitelist file for alevin containing the given dummy cell barcode. 
    """

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Get our model:
    if LibraryModel.has_prebuilt_model(model):
        lb_model = LibraryModel.build_pre_configured_model(model)
    else:
        logger.info(f"Loading model from json file: %s", model)
        lb_model = LibraryModel.from_json_file(model)
    logger.info(f"Using %s: %s", model, lb_model.description)
    model = lb_model

    # Configure process manager:
    # NOTE: We're using processes to overcome the Global Interpreter Lock.
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(threads)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_process_pool = []
    for _ in range(threads):
        p = mp.Process(
            target=_sub_process_work_fn, args=(process_input_data_queue, results, umi_length, model, write_bam)
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
        if force:
            logger.info("Force mode - skipping bam header check for compatibility")
        else:
            # Make sure we're given an input bam file we can work with:
            if not _validate_input_bam(bam_file.header):
                # Bad news - we have to quit.
                # let's try to do it nicely:
                for r in (None,) * threads:
                    process_input_data_queue.put(r)

                # Wait for our input jobs to finish:
                for p in worker_process_pool:
                    p.join()
                sys.exit(1)

        # Get our model:
        if model is None:
            lb_model = LibraryModel.from_json_obj(bam_utils.get_model_from_bam_header(bam_file.header))
        elif model is not None and LibraryModel.has_prebuilt_model(model):
            lb_model = LibraryModel.build_pre_configured_model(model)
        else:
            lb_model = LibraryModel.from_json_file(model)

        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[model])

        # Start output worker:
        res = manager.dict({"num_reads_processed": 0})
        output_worker = mp.Process(
            target=_sub_process_write_fn,
            args=(
                results,
                output_base_name,
                cell_barcode,
                pbar,
                res,
                write_bam,
                out_header
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

    # Write out our CBC whitelist file:
    with open(f"{output_base_name}{__OUT_WHITELIST_FILE_SUFFIX}", "w") as f:
        f.write(f"{cell_barcode}\n")

    logger.info(f"Processed {res['num_reads_processed']} reads.")
    logger.info(f"CBC length: {len(cell_barcode)}.")
    logger.info(f"UMI length: {umi_length}.")
    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)


def _validate_input_bam(input_bam_header):
    """Check that the given input_bam_header contains an `longbow segment` program group."""
    in_bam_header_dict = input_bam_header.to_dict()
    if "PG" not in in_bam_header_dict:
        logger.warning("Could not find PG entry in header.  Cannot confirm that this file is compatible.")
    else:
        found_segment_cmd = False
        for info in [item for item in in_bam_header_dict["PG"]]:
            if "PN" not in info:
                continue
            if info["PN"] == "longbow" and info["ID"].split("-")[1] == "segment":
                found_segment_cmd = True
                break
        if not found_segment_cmd:
            logger.error(
                "Input bam file header does not indicate that it was created by longbow segment.  "
                "This tool requires `longbow segment` reads as input data.")
            return False
    return True


def _get_start_segment_from_list(seg_list, model, read_name):
    """Get the start segment segment from the list of SegmentInfo objects based on the given model.
    If no start segment is found, returns None."""

    # The start segment should be the first matching segment:
    for s in seg_list:
        if s.name in model.start_element_names:
            return s

    logger.warning("Could not process read: %s - No start segment found (start names: %s).",
                   read_name, model.start_element_names)
    return None


def _get_end_segment_from_list(seg_list, model, read_name):
    """Get the end segment segment from the list of SegmentInfo objects based on the given model.
    If no start segment is found, returns None."""

    # The end segment should be the last matching segment, so we
    # iterate from the end to the start of the list:
    for s in reversed(seg_list):
        if s.name in model.end_element_names:
            return s

    logger.warning("Could not process read: %s - No end segment found (end names: %s).",
                   read_name, model.start_element_names)
    return None


def _sub_process_work_fn(in_queue, out_queue, umi_length, array_model, do_bam_out):
    """Function to run in each subprocess.
    Extracts and returns all segments from an input read."""
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
        _, segments = get_segments(read)

        # Get start element position
        # (for MAS-seq it's the 10x adapter)
        start_segment = _get_start_segment_from_list(segments, array_model, read.query_name)
        if start_segment is None:
            continue

        # Get the end element position:
        # (for MAS-seq it's the Poly-a)
        end_segment = _get_end_segment_from_list(segments, array_model, read.query_name)
        if end_segment is None:
            continue

        # Now we grab the bases just after the 10x adapter as the UMI
        # and the bases between the UMI and the poly A for the transcript

        # Note: Positions are inclusive so we must add 1 to the end position to get that base as well:
        umi_start = start_segment.end+1
        umi_end = umi_start + umi_length
        umi_bases = read.query_sequence[umi_start:umi_end]
        umi_quals = "".join([chr(i + 33) for i in read.query_alignment_qualities[umi_start:umi_end]])

        transcript_bases = read.query_sequence[umi_end:end_segment.start]
        transcript_quals = "".join(
            [chr(i + 33) for i in read.query_alignment_qualities[umi_end:end_segment.start]]
        )

        # Place our data on the output queue:
        if do_bam_out:
            out_queue.put(
                tuple([read.query_name, umi_bases, umi_quals, transcript_bases, transcript_quals, read.to_string()])
            )
        else:
            out_queue.put(
                tuple([read.query_name, umi_bases, umi_quals, transcript_bases, transcript_quals])
            )


def _sub_process_write_fn(
    out_queue,
    out_base_name,
    cell_barcode,
    pbar,
    res,
    do_bam_out,
    out_bam_header
):
    """Thread / process fn to write out all our data."""

    try:
        if do_bam_out:
            out_bam_file = pysam.AlignmentFile(f"{out_base_name}.cbc_umi_annotated.bam", "wb", header=out_bam_header)
        with open(f"{out_base_name}{__OUT_READ_FILE_SUFFIX}1.fastq", "w") as mates1_file, \
             open(f"{out_base_name}{__OUT_READ_FILE_SUFFIX}2.fastq", "w") as mates2_file:

            while True:
                # Wait for some output data:
                raw_data = out_queue.get()

                # Check for exit sentinel:
                if raw_data is None:
                    break

                # Unpack data:
                if do_bam_out:
                    read_name, umi_bases, umi_quals, transcript_bases, transcript_quals, read_string = raw_data
                else:
                    read_name, umi_bases, umi_quals, transcript_bases, transcript_quals = raw_data

                # Create mates1 and mates2 records:
                mates_1_record = pysam.FastxRecord(
                    name=read_name,
                    sequence=cell_barcode + umi_bases,
                    quality=(chr(33 + 60) * len(cell_barcode)) + umi_quals
                )
                mates_2_record = pysam.FastxRecord(
                    name=read_name,
                    sequence=transcript_bases,
                    quality=transcript_quals
                )

                # Write out mates1 and mates2 records:
                mates1_file.write(str(mates_1_record))
                mates1_file.write("\n")
                mates2_file.write(str(mates_2_record))
                mates2_file.write("\n")

                if do_bam_out:
                    read = pysam.AlignedSegment.fromstring(
                        read_string, pysam.AlignmentHeader.from_dict(dict())
                    )

                    read.set_tag("CR", cell_barcode)
                    read.set_tag("UR", umi_bases)

                    out_bam_file.write(read)

                # Increment our counters:
                res["num_reads_processed"] += 1
                pbar.update(1)

                # Obligatory log message:
                logger.debug("Processed read: %s", read_name)
    finally:
        if do_bam_out:
            out_bam_file.close()
