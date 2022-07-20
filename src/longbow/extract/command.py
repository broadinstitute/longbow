import logging
import time
import os
import sys

import click
import click_log

import tqdm
import pysam
from construct import *

import longbow.utils.constants
from ..utils import bam_utils
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel

from ..utils.bam_utils import SegmentInfo

from ..annotate.command import get_segments
from ..segment.command import create_simple_delimiters
from ..segment.command import segment_read_with_simple_splitting
from ..segment.command import create_simple_split_array_element

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("extract")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-p",
    "--pbi",
    required=False,
    type=click.Path(exists=True),
    help="BAM .pbi index file",
)
@click.option(
    "-o",
    "--output-bam",
    default="-",
    type=click.Path(exists=False),
    help="extracted bam output.  [default: stdout]",
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
    "--create-barcode-conf-file",
    required=False,
    is_flag=True,
    default=False,
    help=f"Create a barcode confidence score file based on the barcodes in the given model.  "
         f"This only applies for models that have annotation_segments where one such segment "
         f"is annotated into the raw barcode field ({longbow.utils.constants.READ_RAW_BARCODE_TAG})",
)
@click.option(
    "-b",
    "--base-padding",
    default=2,
    required=False,
    show_default=True,
    type=int,
    help="Number of bases to include on either side of the extracted region(s).",
)
@click.option(
    "--leading-adapter",
    default=None,
    required=False,
    type=str,
    help="Adapter preceding the region to extract.  Required if the given model does not name a `coding_region`.",
)
@click.option(
    "--trailing-adapter",
    default=None,
    required=False,
    type=str,
    help="Adapter following the region to extract.  Required if the given model does not name a `coding_region`.",
)
@click.option(
    # CBC + UMI for MAS15 is 26
    "--start-offset",
    default=None,
    required=False,
    show_default=True,
    type=int,
    help="Number of bases to ignore from the extracted region start.  "
         "These bases will not be included in the extracted sequences.  "
         "Required if the given model does not name a `coding_region`.",
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
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, output_bam, force, base_padding, create_barcode_conf_file,
         leading_adapter, trailing_adapter, start_offset, model, input_bam):
    """Extract coding segments from the reads in the given bam.
    The main coding segments are assumed to be labeled as `random` segments or labeled with the same name as the
    `coding_region` in the model (if specified).

    For `random` segment extraction:
    Uses known segments flanking the region to be extracted as markers to indicate
    the start and end of what to extract.

    For `coding_region` segment extraction:
    Looks for any section of the reads labeled with the same name as the `coding_region` in the model, regardless of
    position.
    """

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output file exists:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    logger.info(f"Writing extracted read segments to: {output_bam}")
    logger.info(f"Including {base_padding} flanking bases.")

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("About to Extract segments from %d reads", read_count)

    # Open our input bam file:
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file, \
            tqdm.tqdm(
            desc="Progress",
            unit=" read",
            colour="green",
            file=sys.stderr,
            disable=not sys.stdin.isatty(),
            total=read_count
            ) as pbar:

        # Get our model:
        if model is None:
            lb_model = LibraryModel.from_json_obj(bam_utils.get_model_from_bam_header(bam_file.header))
        elif model is not None and LibraryModel.has_prebuilt_model(model):
            lb_model = LibraryModel.build_pre_configured_model(model)
        else:
            lb_model = LibraryModel.from_json_file(model)

        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        # Validate our command line arguments:
        if lb_model.has_coding_region:
            logger.info(f"Extracting coding region from model {lb_model.name}: {lb_model.coding_region}")
        else:
            # We don't have a model with a coding region.
            # We MUST make sure we have the other inputs that we need:
            missing_args = []

            if leading_adapter is None:
                missing_args.append("--leading-adapter")
            if trailing_adapter is None:
                missing_args.append("--trailing-adapter")
            if start_offset is None:
                missing_args.append("--start-offset")

            if len(missing_args) > 0:
                message = f"ERROR: the following arguments must be specified when using a model that has no " \
                        f"`coding_region`: {', '.join(missing_args)}"
                logger.critical(message)
                raise RuntimeError(message)

            logger.info(f"Extracting `{longbow.utils.constants.RANDOM_SEGMENT_NAME}` segments between {leading_adapter} and {trailing_adapter}.")
            logger.info(f"Ignoring the first {start_offset} bases from extracted read segments.")

        # Set up our barcode confidence file here:
        barcode_conf_file = None
        if create_barcode_conf_file:
            if lb_model.has_cell_barcode_annotation:
                logger.info(f"Creating barcode confidence file: {longbow.utils.constants.BARCODE_CONF_FILE_NAME}")
                barcode_conf_file = open(longbow.utils.constants.BARCODE_CONF_FILE_NAME, 'w')
            else:
                logger.warning(f"Model does not have a barcode output, but barcode creation flag was given.  "
                               f"Barcode confidence file will NOT be created.")

        # Get our delimiters here just in case we have to split the reads later:
        delimiters = create_simple_delimiters(lb_model)

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header)

        # Setup output files:
        with pysam.AlignmentFile(output_bam, "wb", header=out_header) as extracted_bam_file:

            num_reads = 0
            num_reads_with_extracted_segments = 0
            num_segments_extracted = 0
            num_segments_skipped = 0

            for read in bam_file:

                # Get our read segments:
                try:
                    _, segments = get_segments(read)
                except KeyError:
                    logger.error(f"Input bam file does not contain longbow segmented reads!  "
                                 f"No {longbow.utils.constants.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                    sys.exit(1)

                # Check if the read is already segmented:
                if (not read.has_tag(longbow.utils.constants.READ_IS_SEGMENTED_TAG)) or (not read.get_tag(
                        longbow.utils.constants.READ_IS_SEGMENTED_TAG)):
                    # The read is not segmented.  We should segment it first and then go through our segments one by one:
                    logger.debug(f"Read must be segmented prior to extraction: {read.query_name}")
                    segmented_reads = []
                    segment_bounds_tuples = segment_read_with_simple_splitting(read, delimiters, segments)
                    for prev_delim_name, delim_name, start_coord, end_coord in segment_bounds_tuples:
                        segmented_reads.append(
                            create_simple_split_array_element(
                                delim_name, end_coord, lb_model, prev_delim_name, read, segments, start_coord
                            )
                        )
                else:
                    # Read is already segmented, so we don't have to do anything else:
                    segmented_reads = [read]

                extracted_segment = False
                for r in segmented_reads:
                    segments = _get_segments_only(r)
                    # If our model has a coding region, we should use it:
                    if lb_model.has_coding_region:
                        es, nsx, nss = \
                            _extract_region_from_model_with_coding_region(r, segments, lb_model, base_padding, extracted_bam_file, barcode_conf_file)
                    else:
                        # We do not have a coding section, so we'll need to do things the old fashioned way:
                        es, nsx, nss = \
                            _extract_region_from_model_without_coding_region(r, segments, leading_adapter, trailing_adapter, start_offset, base_padding, extracted_bam_file, barcode_conf_file)

                    # Track our stats:
                    extracted_segment |= es
                    num_segments_extracted += nsx
                    num_segments_skipped += nss

                num_reads += 1
                if extracted_segment:
                    num_reads_with_extracted_segments += 1

                pbar.update(1)

    # Close our barcode file if it was opened:
    if barcode_conf_file is not None:
        barcode_conf_file.close()

    # Calc some stats:
    pct_reads_with_extracted_segments = 100 * num_reads_with_extracted_segments / num_reads if num_reads > 0 else 0
    segs_per_read = num_segments_extracted / num_reads if num_reads > 0 else 0

    # Yell at the user:
    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)
    logger.info(f"Total # Reads Processed: %d", num_reads)
    logger.info(f"# Reads Containing Extracted Segments: %d (%2.2f%%)",
                num_reads_with_extracted_segments, pct_reads_with_extracted_segments)
    logger.info(f"Total # Segments Extracted: %d", num_segments_extracted)
    logger.info(f"Total # Segments Skipped: %d", num_segments_skipped)
    logger.info(f"# Segments extracted per read: %2.2f", segs_per_read)


def _get_segments_only(read):
    """Get the segments corresponding to a particular read by reading the segments tag information."""
    return [
        SegmentInfo.from_tag(s) for s in read.get_tag(longbow.utils.constants.SEGMENTS_TAG).split(
            longbow.utils.constants.SEGMENT_TAG_DELIMITER)
    ]


def _extract_region_from_model_without_coding_region(array_element_read, segments, leading_adapter,
                                                     trailing_adapter, start_offset, base_padding, extracted_bam_file,
                                                     barcode_conf_file):
    # Get our marker segments:
    start_marker_list = [(i, s) for i, s in enumerate(segments) if s.name == leading_adapter]
    end_marker_list = [(i, s) for i, s in enumerate(segments) if s.name == trailing_adapter]

    if len(start_marker_list) != len(end_marker_list):
        logger.warning(f"Found %d start markers and %d end markers.  Only looking at first %d pairs.  "
                       f"(starts: %s, ends: %s)",
                       len(start_marker_list), len(end_marker_list),
                       min(len(start_marker_list), len(end_marker_list)),
                       " ".join([f"{i}:{s.name}" for i, s in start_marker_list]),
                       " ".join([f"{i}:{s.name}" for i, s in end_marker_list]))

    extracted_segment = False

    num_segments_skipped = 0
    num_segments_extracted = 0

    # Go through each marker pair and do the extraction:
    for s_info, e_info in zip(start_marker_list, end_marker_list):

        si, start_marker = s_info
        ei, end_marker = e_info

        # Does the start marker come before the end marker and do we have exactly one random segment in
        # between them?
        if (start_marker.end < end_marker.start) and (ei - si == 2) and (
                segments[si + 1].name == "random"):
            # We have a valid segment to extract:
            logger.debug("Found a segment to extract: %s: %s", array_element_read.query_name, segments[si + 1])

            # Create an AlignedSegment to output:
            aligned_segment = _create_extracted_aligned_segment(
                array_element_read, segments[si + 1], start_offset, base_padding, barcode_conf_file
            )

            if aligned_segment:
                extracted_bam_file.write(aligned_segment)
                num_segments_extracted += 1
                extracted_segment = True
            else:
                num_segments_skipped += 1
        else:
            if start_marker.end >= end_marker.start:
                logger.warning(
                    "Read %s: start marker segment (i=%d) occurs at or after end segment (i=%d):"
                    " %d >= %d.  Skipping segment.",
                    array_element_read.query_name, si, ei, start_marker.end, end_marker.start)
            elif ei - si != 2:
                logger.warning(
                    "Read %s: start segment (i=%d) and end segment (i=%d) have more than one "
                    "segment between them.  Skipping segment.", array_element_read.query_name, si, ei)
            elif segments[si + 1].name != "random":
                logger.warning("Read %s: segment between start segment (i=%d) and end segment (i=%d) "
                               "is not a random segment.  Skipping segment.", array_element_read.query_name, si, ei)
            num_segments_skipped += 1

    return extracted_segment, num_segments_extracted, num_segments_skipped


def _extract_region_from_model_with_coding_region(array_element_read, segments, m, base_padding, extracted_bam_file, barcode_conf_file):
    # Get our marker segments:
    # NOTE: Even though this is a list, we only expect one segment to be here:
    extraction_segments = [s for s in segments if s.name == m.coding_region]

    if len(extraction_segments) == 0:
        logger.debug(f"Did not find coding region in read: %s: %s", array_element_read.query_name, segments)
        return False, 0, 1

    else:
        if len(extraction_segments) > 1:
            logger.debug(f"Found %d coding regions in read %s: %s.  Only looking at first coding region.  All Segments: %s",
                         len(extraction_segments), array_element_read.query_name, extraction_segments, segments)

        # Create an AlignedSegment to output:
        aligned_segment = _create_extracted_aligned_segment(
            array_element_read, extraction_segments[0], 0, base_padding, barcode_conf_file
        )

        if aligned_segment:
            extracted_bam_file.write(aligned_segment)
            return True, 1, 0
        else:
            return False, 0, 1


def _create_extracted_aligned_segment(read, seg_to_extract, start_offset, base_padding, barcode_conf_file):
    """Create a pysam.AlignedSegment object to store the information from the extracted bases."""

    start_coord = seg_to_extract.start + start_offset - base_padding
    end_coord = seg_to_extract.end + base_padding

    # Bounds check our coords:
    if start_coord < 0:
        logger.debug("Calculated start for %s would start before read begins.  Setting to 0.", read.query_name)
        start_coord = 0
    elif start_coord >= len(read.query_sequence):
        logger.warning("Start coord for %s would start after read.  Cannot process.", read.query_name)
        return None

    if end_coord < 0:
        logger.warning("End coord for %s would start before read.  Cannot process.", read.query_name)
        return None
    elif end_coord >= len(read.query_sequence):
        logger.debug("Calculated end for %s would start after read ends.  Setting to 0.", read.query_name)
        end_coord = len(read.query_sequence)-1

    if end_coord <= start_coord:
        logger.warning("Start coord for %s would start at or after end coord.  Cannot process.", read.query_name)
        return None

    # Create our segment:
    a = pysam.AlignedSegment()
    a.query_name = read.query_name
    # Add one to end_coord because coordinates are inclusive:
    a.query_sequence = f"{read.query_sequence[start_coord:end_coord+1]}"
    a.query_qualities = read.query_alignment_qualities[start_coord: end_coord + 1]
    a.tags = read.get_tags()
    a.flag = 4  # unmapped flag
    a.mapping_quality = 255
    a.set_tag(longbow.utils.constants.READ_ALTERED_NAME_TAG, f"{read.query_name}/{start_coord}_{end_coord}")

    # Write our barcode confidence to the file if we have to:
    if barcode_conf_file is not None and a.has_tag(longbow.utils.constants.READ_BARCODE_CONF_FACTOR_TAG):
        barcode_conf_file.write(f"{a.get_tag(longbow.utils.constants.READ_RAW_BARCODE_TAG)}\t{a.get_tag(longbow.utils.constants.READ_BARCODE_CONF_FACTOR_TAG)}\n")

    return a
