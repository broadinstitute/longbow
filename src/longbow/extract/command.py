import logging
import time
import os
import sys

import click
import click_log

import tqdm
import pysam
from construct import *

from ..utils import bam_utils
from ..annotate.command import get_segments

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
    "--out-file",
    default=".",
    required=True,
    type=str,
    help="Output file name.",
)
@click.option(
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
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
    default="10x_Adapter",
    required=False,
    show_default=True,
    type=str,
    help="Adapter preceding the region to extract.",
)
@click.option(
    "--trailing-adapter",
    default="Poly_A",
    required=False,
    show_default=True,
    type=str,
    help="Adapter following the region to extract.",
)
@click.option(
    "--start-offset",
    default=16+10,  # CBC + UMI for MAS15 is the default
    required=False,
    show_default=True,
    type=int,
    help="Number of bases to ignore from the extracted region start.  "
         "These bases will not be included in the extracted sequences.",
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, out_file, force, base_padding, leading_adapter, trailing_adapter, start_offset, input_bam):
    """Extract coding segments from the reads in the given bam.
    The main coding segments are assumed to be labeled as `random` segments.
    Uses known segments flanking the region to be extracted as markers to indicate
    the start and end of what to extract."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("About to Extract segments from %d reads", read_count)

    # Check to see if the output file exists:
    bam_utils.check_for_preexisting_files(out_file, exist_ok=force)

    # TODO: We don't need to check our model right now because our models are SO similar.
    #       This will need to be fixed when we start using more exotic models.

    logger.info(f"Writing extracted read segments to: {out_file}")
    logger.info(f"Extracting `random` segments between {leading_adapter} and {trailing_adapter}.")
    logger.info(f"Ignoring the first {start_offset} bases from extracted read segments.")
    logger.info(f"Including {base_padding} flanking bases.")

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

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header)

        # Setup output files:
        with pysam.AlignmentFile(out_file, "wb", header=out_header) as extracted_bam_file:

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
                                 f"No {bam_utils.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                    sys.exit(1)

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

                # Go through each marker pair and do the extraction:
                for s_info, e_info in zip(start_marker_list, end_marker_list):

                    si, start_marker = s_info
                    ei, end_marker = e_info

                    # Does the start marker come before the end marker and do we have exactly one random segment in
                    # between them?
                    if (start_marker.end < end_marker.start) and (ei - si == 2) and (segments[si + 1].name == "random"):
                        # We have a valid segment to extract:
                        logger.debug("Found a segment to extract: %s: %s", read.query_name, segments[si + 1])

                        # Create an AlignedSegment to output:
                        aligned_segment = _create_extracted_aligned_segment(
                            read, segments[si + 1], start_offset, base_padding
                        )

                        if aligned_segment:
                            extracted_bam_file.write(aligned_segment)
                            num_segments_extracted += 1
                            extracted_segment = True
                        else:
                            num_segments_skipped += 1
                    else:
                        if start_marker.end >= end_marker.start:
                            logger.warning("Read %s: start marker segment (i=%d) occurs at or after end segment (i=%d):"
                                           " %d >= %d.  Skipping segment.",
                                           read.query_name, si, ei, start_marker.end, end_marker.start)
                        elif ei - si != 2:
                            logger.warning("Read %s: start segment (i=%d) and end segment (i=%d) have more than one "
                                           "segment between them.  Skipping segment.", read.query_name, si, ei)
                        elif segments[si + 1].name != "random":
                            logger.warning("Read %s: segment between start segment (i=%d) and end segment (i=%d) "
                                           "is not a random segment.  Skipping segment.", read.query_name, si, ei)
                        num_segments_skipped += 1

                pbar.update(1)
                num_reads += 1
                if extracted_segment:
                    num_reads_with_extracted_segments += 1

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


def _create_extracted_aligned_segment(read, seg_to_extract, start_offset, base_padding):
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
    a.query_name = (
        f"{read.query_name}/{start_coord}_{end_coord}"
    )
    # Add one to end_coord because coordinates are inclusive:
    a.query_sequence = f"{read.query_sequence[start_coord:end_coord+1]}"
    a.query_qualities = read.query_alignment_qualities[start_coord: end_coord + 1]
    a.tags = read.get_tags()
    a.flag = 4  # unmapped flag
    a.mapping_quality = 255

    return a