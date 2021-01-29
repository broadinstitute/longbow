import logging
import sys
import math

from collections import namedtuple

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp
import concurrent.futures

from ..utils.model import build_default_model
from ..utils.model import array_element_structure
from ..utils.model import reverse_complement
from ..utils.model import annotate


# Named tuple to store alignment information:
class SegmentInfo(namedtuple("SegmentInfo", ["name", "start", "end"])):
    def __str__(self):
        return f"SegmentInfo({self.to_tag()})"

    def to_tag(self):
        return f"{self.name}:{self.start}-{self.end}"


logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command(name="segment")
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
@click.option(
    "-p",
    "--split-bam",
    required=False,
    type=click.Path(exists=False),
    help="segmented bam output",
)
@click.option(
    "-s",
    "--do-simple-splitting",
    required=False,
    is_flag=True,
    default=False,
    help="Do splitting of reads based on splitter delimiters, rather than whole array structure."
    "This splitting will cause delimiter sequences to be repeated in each read they bound.",
)
@click.argument("input-bam", type=click.Path(exists=True))
def main(model, threads, output_bam, split_bam, do_simple_splitting, input_bam):
    """Apply annotation and segmentation model to BAM file"""
    logger.info("annmas: segment started")

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} thread(s)")

    m = build_default_model()
    if model is not None:
        m.from_yaml(model)
        logger.info(f"Using pretrained annotation model {model}")
    else:
        logger.info("Using default annotation model")

    if do_simple_splitting:
        logger.info("Using simple splitting mode.")
    else:
        logger.info("Using bounded region splitting mode.")

    num_reads_annotated = 0
    num_sections = 0
    num_segments = 0

    # Placeholder for split out file:
    # NOTE: There's probably a more pythonic way to do this...
    split_bam_out = None

    try:
        pysam.set_verbosity(0)  # silence message about the .bai file not being found
        with pysam.AlignmentFile(
            input_bam, "rb", check_sq=False, require_index=False
        ) as bam_file:

            with concurrent.futures.ThreadPoolExecutor(
                max_workers=threads
            ) as executor, tqdm.tqdm(
                desc="Progress", unit=" read", colour="green", file=sys.stdout
            ) as pbar:

                future_to_segmented_read = {
                    executor.submit(_segment_read, r, m): r for r in bam_file
                }

                # Get our header from the input bam file:
                out_bam_header_dict = bam_file.header.to_dict()

                # Add our program group to it:
                pg_dict = {
                    "ID": "annmas-segment-0.0.1",
                    "PN": "annmas",
                    "VN": "0.0.1",
                    "DS": "Apply annotation and segmentation model to BAM file.",
                    "CL": " ".join(sys.argv),
                }
                if "PG" in out_bam_header_dict:
                    out_bam_header_dict["PG"].append(pg_dict)
                else:
                    out_bam_header_dict["PG"] = [pg_dict]

                if split_bam:
                    split_bam_out = pysam.AlignmentFile(
                        split_bam, "wb", header=out_bam_header_dict
                    )

                with pysam.AlignmentFile(
                    output_bam, "wb", header=out_bam_header_dict
                ) as out_bam_file:

                    for future in concurrent.futures.as_completed(
                        future_to_segmented_read
                    ):
                        read = future_to_segmented_read[future]
                        try:
                            path, logp = future.result()
                        except Exception as ex:
                            logger.error("%r generated an exception: %s", read, ex)
                        else:
                            # Condense the output annotations so we can write them out with indices:
                            segments = _collapse_annotations(path)

                            # Obligatory log message:
                            logger.debug(
                                "Path for read %s (%2.2f): %s",
                                read.query_name,
                                logp,
                                segments,
                            )

                            # Set our tag and write out the read to the annotated file:
                            read.set_tag("SG", "|".join([s.to_tag() for s in segments]))
                            out_bam_file.write(read)

                            # Write out our segmented reads if we want them:
                            if split_bam_out:
                                num_segments += _write_segmented_read(
                                    read, segments, do_simple_splitting, split_bam_out
                                )
                            # Increment our counters:
                            num_reads_annotated += 1
                            num_sections += len(segments)

                            pbar.update(1)

        logger.info("annmas: segment finished.")
        logger.info(
            f"annmas: segmented {num_reads_annotated} reads into {num_sections} sections."
        )

    finally:
        if split_bam_out:
            split_bam_out.close()
            logger.info(f"annmas: wrote {num_segments} segments.")


def _write_segmented_read(read, segments, do_simple_splitting, bam_out):
    """Split and write out the segments of each read to the given bam output file
    :param read: A pysam.AlignedSegment object containing a read that has been segmented.
    :param segments: A list of SegmentInfo objects representing the segments of the given reads.
    :param do_simple_splitting: Flag to control how reads should be split.
                                If True, will use simple delimiters.
                                If False, will require reads to appear as expected in model.array_element_structure.
    :param bam_out: An open pysam.AlignmentFile ready to write out data.
    :return: the number of segments written.
    """

    if do_simple_splitting:
        # Create the sections on which we want to split.
        num_required_delimiters = 2

        delimiters = list()
        delimiters.append(tuple(array_element_structure[0][-num_required_delimiters:]))
        for i, structure in enumerate(array_element_structure[1:], start=1):
            # If it's the second element then we append the delmiters to those from the first.
            # This is simple splitting, after all.
            if i == 1:
                delimiters[0] = delimiters[0] + structure[0:num_required_delimiters]
            else:
                delimiters.append(tuple(structure[0:num_required_delimiters]))

        # Now we have our delimiter list.
        # We need to go through our segments and split them up.
        # Note: we assume each delimiter can occur only once.
        delimiter_match_matrix = [0 for _ in delimiters]
        delimiter_start_segments = [None for _ in delimiters]

        # We have to store the end segments as tuples in case we want to use more than one
        # segment as a delimiter:
        delimiter_end_segment_tuples = [None for _ in delimiters]

        # We do it this way so we iterate over the segments once and the delimiters many times
        # under the assumption the segment list is much longer than the delimiters.
        for seg in segments:
            # at each position go through our delimiters and track whether we're a match:
            for i, dmi in enumerate(delimiter_match_matrix):
                try:
                    if seg.name == delimiters[i][dmi]:
                        if delimiter_match_matrix[i] == 0:
                            delimiter_start_segments[i] = seg
                        delimiter_match_matrix[i] += 1
                        if delimiter_match_matrix[i] == len(delimiters[i]):
                            # We've got a full and complete match!
                            # We store the end segment we found.
                            delimiter_end_segment_tuples[i] = (
                                delimiter_start_segments[i],
                                seg,
                            )
                    else:
                        # No match at the position we expected!
                        # We need to reset our count and start segment:
                        delimiter_match_matrix[i] = 0
                        delimiter_start_segments[i] = None
                        delimiter_end_segment_tuples[i] = None
                except IndexError:
                    # We're out of range of the delimiter.
                    # This means we've hit the end and actually have a match, so we can ignore this error.
                    pass

        # OK, we've got a handle on our delimiters, so now we just need to split the sequence.
        # We do so by looking at which end delimiters are filled in, and grabbing their complementary start delimiter:

        # Sort and filter our delimiters:
        seg_delimiters = [
            (delim_index, start, end_tuple)
            for delim_index, (start, end_tuple) in enumerate(
                zip(delimiter_start_segments, delimiter_end_segment_tuples)
            )
            if end_tuple is not None
        ]
        seg_delimiters.sort(key=lambda dtuple: dtuple[1].start)

        cur_read_base_index = 0
        prev_delim_name = "START"

        for i, (di, start_seg, end_seg_tuple) in enumerate(seg_delimiters):

            start_coord = cur_read_base_index
            end_coord = end_seg_tuple[1].end
            delim_name = "/".join(delimiters[di])

            # Write our segment here:
            _write_split_array_element(
                bam_out,
                start_coord,
                end_coord,
                read,
                segments,
                delim_name,
                prev_delim_name,
            )

            cur_read_base_index = end_seg_tuple[0].start
            prev_delim_name = delim_name

        # Now we have to write out the last segment:
        start_coord = cur_read_base_index
        end_coord = len(read.query_sequence)
        delim_name = "END"

        _write_split_array_element(
            bam_out, start_coord, end_coord, read, segments, delim_name, prev_delim_name
        )

        return len(seg_delimiters)
    else:
        # Here we're doing bounded region splitting.
        # This requires each read to conform to the expected read structure as defined in the model.
        # The process is similar to what is done above for simple splitting.

        # Create our delimiter list.
        delimiters = array_element_structure

        # We need to go through our segments and split them up.
        # Note: we assume each full array element can occur only once and they do not overlap.
        delimiter_match_matrix = [0 for _ in delimiters]
        delimiter_segments = [list() for _ in delimiters]
        delimiter_found = [False for _ in delimiters]
        delimiter_score = [0 for _ in delimiters]

        # We define some scoring increments here:
        match_val = 2
        indel_val = 1

        # We do it this way so we iterate over the segments once and the delimiters many times
        # under the assumption the segment list is much longer than the delimiters.
        for seg in segments:
            # at each position go through our delimiters and track whether we're a match:
            for i, dmi in enumerate(delimiter_match_matrix):
                if not delimiter_found[i]:
                    try:
                        if seg.name == delimiters[i][dmi]:
                            delimiter_match_matrix[i] += 1
                            delimiter_segments[i].append(seg)
                            delimiter_score[i] += match_val
                        else:
                            found = False
                            # Only look ahead if we already have a partial match:
                            if delimiter_match_matrix[i] != 0:
                                # Here we "peek ahead" so we an look at other delimiters after the current one just in
                                # case we're missing a delimiter / segment.  This will impact the "score" but allow for
                                # fuzzy matching.
                                for peek_ahead in range(
                                    1, len(delimiters[i]) - delimiter_match_matrix[i]
                                ):
                                    if seg.name == delimiters[i][dmi + peek_ahead]:
                                        delimiter_match_matrix[i] += 1 + peek_ahead
                                        delimiter_segments[i].append(seg)
                                        delimiter_score[i] += indel_val
                                        found = True
                                        break

                            if not found:
                                # No match at the position we expected!
                                # We need to reset our count and start segment:
                                delimiter_match_matrix[i] = 0
                                delimiter_segments[i] = list()
                                delimiter_score[i] = 0

                        # Make sure we mark ourselves as done if we're out of info:
                        if delimiter_match_matrix[i] == len(delimiters[i]):
                            delimiter_found[i] = True

                    except IndexError:
                        # We're out of range of the delimiter.
                        # This means we've hit the end and actually have a match.
                        delimiter_found[i] = True

        # Now we have our segments as described by our model.
        # We assume they don't overlap and we write them out:
        for i, seg_list in enumerate(delimiter_segments):
            if delimiter_found[i]:

                start_seg = seg_list[0]
                end_seg = seg_list[-1]

                start_coord = start_seg.start
                end_coord = end_seg.end
                start_delim_name = seg_list[0].name
                end_delim_name = seg_list[-1].name

                # Write our segment here:
                _write_split_array_element(
                    bam_out,
                    start_coord,
                    end_coord,
                    read,
                    seg_list,
                    end_delim_name,
                    start_delim_name,
                )

        # Return the number of array elements.
        # NOTE: this works because booleans are a subset of integers in python.
        return sum(delimiter_found)


def _write_split_array_element(
    bam_out, start_coord, end_coord, read, segments, delim_name, prev_delim_name
):
    """Write out an individual array element that has been split out according to the given coordinates."""
    a = pysam.AlignedSegment()
    a.query_name = (
        f"{read.query_name}_{start_coord}-{end_coord}_{prev_delim_name}-{delim_name}"
    )
    a.query_sequence = f"{read.query_sequence[start_coord:end_coord]}"
    a.query_qualities = read.query_alignment_qualities[start_coord:end_coord]
    a.tags = read.get_tags()
    a.flag = 4  # unmapped flag
    a.mapping_quality = 255

    # Set our segments tag to only include the segments in this read:
    a.set_tag(
        "SG",
        ",".join([s.to_tag() for s in segments if start_coord <= s.start <= end_coord]),
    )

    bam_out.write(a)


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


def _segment_read(read, um):
    flogp = -math.inf
    fppath = []
    for seq in [read.query_sequence, reverse_complement(read.query_sequence)]:
        logp, ppath = annotate(um, seq)

        if logp > flogp:
            flogp = logp
            fppath = ppath

    return fppath, flogp
