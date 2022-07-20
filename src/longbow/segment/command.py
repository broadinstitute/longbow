import logging
import sys
import itertools
import re

import time

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp
import numpy as np

import longbow.utils.constants
from ..utils import bam_utils
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel

from ..annotate.command import SegmentInfo
from ..annotate.command import get_segments


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("segment")
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
    help="segment-annotated bam output  [default: stdout]",
)
@click.option(
    "-b",
    "--create-barcode-conf-file",
    required=False,
    is_flag=True,
    default=False,
    help=f"Create a barcode confidence score file based on the barcodes in the given model.  "
         f"This only applies for models that have annotation_segments where one such segment "
         f"is annotated into the raw barcode field ({longbow.utils.constants.READ_RAW_BARCODE_TAG})",
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
    '-i',
    '--ignore-cbc-and-umi',
    is_flag=True,
    default=False,
    show_default=True,
    help="Do not require passing reads to have CBC and UMI."
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(threads, output_bam, create_barcode_conf_file, model, ignore_cbc_and_umi, force, input_bam):
    """Segment pre-annotated reads from an input BAM file."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    logger.info("Using simple splitting mode.")

    # Configure process manager:
    # NOTE: We're using processes to overcome the Global Interpreter Lock.
    manager = mp.Manager()
    process_input_data_queue = manager.Queue(threads)
    results = manager.Queue()

    # Start worker sub-processes:
    worker_process_pool = []
    for _ in range(threads):
        p = mp.Process(
            target=_sub_process_work_fn, args=(process_input_data_queue, results)
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

        if ignore_cbc_and_umi:
            logger.info("Ignoring CBC / UMI - all split elements will be written.")
        else:
            # Check to see if the model has a CBC / UMI.  If not, we should turn on this flag
            # and warn the user:

            has_cbc_or_umi_annotation = False
            if lb_model.has_umi_annotation:
                logger.info("Model has UMI annotation.")
            else:
                has_cbc_or_umi_annotation = True

            if lb_model.has_cell_barcode_annotation:
                logger.info("Model has Cell Barcode annotation.")
            else:
                has_cbc_or_umi_annotation = True

            if not has_cbc_or_umi_annotation:
                logger.warning("Model does not have Cell Barcode or UMI tags.  All segments will be emitted.")
                ignore_cbc_and_umi = True

        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[lb_model])

        # Start output worker:
        res = manager.dict({"num_reads_segmented": 0, "num_segments": 0})
        output_worker = mp.Process(
            target=_sub_process_write_fn,
            args=(
                results,
                out_header,
                output_bam,
                pbar,
                create_barcode_conf_file,
                lb_model,
                res,
                ignore_cbc_and_umi
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
        f"Segmented {res['num_reads_segmented']} reads with {res['num_segments']} total segments."
    )
    num_reads = res['num_reads_segmented']
    num_segmented = res['num_segments']
    logger.info(f"MAS-seq gain factor: {num_segmented/num_reads:.02f}x")
    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)


def _sub_process_work_fn(in_queue, out_queue):
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

        # Process and place our data on the output queue:
        out_queue.put(get_segments(read))


def _sub_process_write_fn(
    out_queue,
    out_bam_header,
    out_bam_file_name,
    pbar,
    create_barcode_conf_file,
    model,
    res,
    ignore_cbc_and_umi
):
    """Thread / process fn to write out all our data."""

    delimiters = create_simple_delimiters(model)
    logger.debug(f"Delimiter sequences for simple delimiters: %s", delimiters)

    barcode_conf_file = None
    if create_barcode_conf_file:
        if model.has_cell_barcode_annotation:
            logger.info(f"Creating barcode confidence file: {longbow.utils.constants.BARCODE_CONF_FILE_NAME}")
            barcode_conf_file = open(longbow.utils.constants.BARCODE_CONF_FILE_NAME, 'w')
        else:
            logger.warning(f"Model does not have a barcode output, but barcode creation flag was given.  "
                           f"Barcode confidence file will NOT be created.")

    model_annotates_cbc = model.has_cell_barcode_annotation
    model_annotates_umi = model.has_umi_annotation

    with pysam.AlignmentFile(
        out_bam_file_name, "wb", header=out_bam_header
    ) as out_bam_file:

        while True:
            # Wait for some output data:
            raw_data = out_queue.get()

            # Check for exit sentinel:
            if raw_data is None:
                break

            # Unpack data:
            read, segments = raw_data
            read = pysam.AlignedSegment.fromstring(read, out_bam_header)

            # Obligatory log message:
            logger.debug(
                "Segments for read %s: %s",
                read.query_name,
                segments,
            )

            # Write out our segmented reads:
            res["num_segments"] += _write_segmented_read(
                model,
                read,
                segments,
                delimiters,
                out_bam_file,
                barcode_conf_file,
                ignore_cbc_and_umi,
                model_annotates_cbc,
                model_annotates_umi
            )

            # Increment our counters:
            res["num_reads_segmented"] += 1
            pbar.update(1)

    if barcode_conf_file is not None:
        barcode_conf_file.close()


def create_simple_delimiters(model, num_seqs_from_each_array_element=1):
    """Create delimiters for simple splitting.

    :param model: Model to use for array structure.
    :param num_seqs_from_each_array_element: The number of sequences from each array element to use as delimiters.
    The more sequences used, the more specific the delimiter sequences will be.  A value of 1 will result in delimiter
    tuples that are each 2 sequences long - the last sequence from the leading array element and the first from the
    trailing array element.

    :return: A list of delimiter tuples, each tuple contains names of sequences from our model in the order in which
             they appear in the sequences.
    """

    # We keep a list of delimiter tuples that will serve as the breakpoints in our sequences:
    delimiters = list()

    # Collect the delimiters by grabbing the inner bounding regions of each array element in the order
    # that we expect to see them in the overall library:
    for i in range(1, len(model.array_element_structure)):
        delimiters.append(
            tuple(model.array_element_structure[i - 1][-num_seqs_from_each_array_element:])
            + tuple(model.array_element_structure[i][0:num_seqs_from_each_array_element])
        )

    return delimiters


def segment_read_with_simple_splitting(read, delimiters, segments=None):
    """Segments the given read using the simple splitting algorithm.

    NOTE: Assumes that all given data are in the forward direction.

    :param read: A pysam.AlignedSegment object containing a read that has been segmented.
    :param delimiters: A list of tuples containing the names of delimiter sequences to use to split the given read.
    :param segments: None or A list of SegmentInfo objects representing the segments of the given reads.
                     If None, will be populated by getting segments from the given read.
    :return: a list of tuples containing the split segments in the given read
    """

    if not segments:
        segments = get_segments(read)

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

    # Now we can go through our split sequences and write them out.
    # However we need to account for the "extra" sequences leading and trailing each array element.
    # These are caused by the delimiter sequences taking the first and second halves from different array elements
    # to increase splitting sensitivity.

    cur_read_base_index = 0
    prev_delim_name = "START"

    segment_tuples = []

    for i, (di, start_seg, end_seg_tuple) in enumerate(seg_delimiters):
        start_coord = cur_read_base_index
        end_coord = end_seg_tuple[0].end
        delim_name = end_seg_tuple[0].name

        segment_tuples.append(
            tuple([prev_delim_name, delim_name, start_coord, end_coord])
        )

        cur_read_base_index = end_seg_tuple[1].start
        prev_delim_name = end_seg_tuple[1].name

    # Now we have to write out the last segment:
    seg_start_coord = cur_read_base_index
    # Subtract 1 for 0-based inclusive coords:
    seg_end_coord = len(read.query_sequence) - 1

    start_coord = seg_start_coord
    end_coord = seg_end_coord

    delim_name = "END"

    segment_tuples.append(
        tuple([prev_delim_name, delim_name, start_coord, end_coord])
    )

    return segment_tuples


def segment_read_with_bounded_region_algorithm(read, model, segments=None):
    """Segments the given read using the simple splitting algorithm.

    NOTE: Assumes that all given data are in the forward direction.

    :param read: A pysam.AlignedSegment object containing a read that has been segmented.
    :param model: The model to use for the array segment information.
    :param segments: None or A list of SegmentInfo objects representing the segments of the given reads.
                     If None, will be populated by getting segments from the given read.
    :return: Two tuples:
             delimiter_found: a tuple of booleans indicating if a delimiter in `delimiter_segments` is the start of a
                              segment.
             delimiter_segments: a tuple of lists of SegmentInfo representing the delimiters for each segment in the
                                 read.
    """

    # Get our segments if we have to:
    if not segments:
        segments = get_segments(read)

    # Create our delimiter list:
    delimiters = model.array_element_structure

    # Modify the delimiter sequences for this library.
    # NOTE: this is an artifact of some library preparation steps and should be removed in the next
    #       iteration of the library.
    # We must change to a list here so we cna modify it:
    delimiters = list(delimiters)
    # The first sequence doesn't actually have an 'A' element:
    delimiters[0] = delimiters[0][1:]
    # The last sequence doesn't actually have a 'P' element:
    delimiters[-1] = delimiters[-1][:-1]

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

    return tuple(delimiter_found), tuple(delimiter_segments)


def _write_segmented_read(
    model, read, segments, delimiters, bam_out, barcode_conf_file, ignore_cbc_and_umi, model_annotates_cbc, model_annotates_umi
):
    """Split and write out the segments of each read to the given bam output file.

    NOTE: Assumes that all given data are in the forward direction.

    :param model: The model to use for the array segment information.
    :param read: A pysam.AlignedSegment object containing a read that has been segmented.
    :param segments: A list of SegmentInfo objects representing the segments of the given reads.
    :param delimiters: A list of tuples containing the names of delimiter sequences to use to split the given read.
    :param bam_out: An open pysam.AlignmentFile ready to write out data.
    :param barcode_conf_file: An open file ready to write out the barcodes and confidence scores.
    :param ignore_cbc_and_umi: Boolean indicating whether we should ignore CBC/UMI tags.
    :return: the number of segments written.
    """

    segment_bounds_tuples = segment_read_with_simple_splitting(read, delimiters, segments)

    sri = 1
    for prev_delim_name, delim_name, start_coord, end_coord in segment_bounds_tuples:
        # Write our segment here:
        _write_split_array_element(
            model,
            bam_out,
            barcode_conf_file,
            start_coord,
            end_coord,
            read,
            segments,
            delim_name,
            prev_delim_name,
            ignore_cbc_and_umi,
            model_annotates_cbc,
            model_annotates_umi,
            sri
        )

        sri += 1

    return len(segment_bounds_tuples)


def _transform_to_rc_coords(start, end, read_length):
    """Transforms the given start and end coordinates into the RC equivalents using the given read_length."""
    return read_length - end - 1, read_length - start - 1


def _write_split_array_element(
    model,
    bam_out,
    barcode_conf_file,
    start_coord,
    end_coord,
    read,
    segments,
    delim_name,
    prev_delim_name,
    ignore_cbc_and_umi,
    model_annotates_cbc,
    model_annotates_umi,
    split_read_index
):
    """Write out an individual array element that has been split out according to the given coordinates."""
    a = create_simple_split_array_element(delim_name, end_coord, model, prev_delim_name, read, segments, start_coord, split_read_index)

    # Write our barcode confidence to the file if we have to:
    if barcode_conf_file is not None and a.has_tag(longbow.utils.constants.READ_BARCODE_CONF_FACTOR_TAG):
        barcode_conf_file.write(f"{a.get_tag(longbow.utils.constants.READ_RAW_BARCODE_TAG)}\t{a.get_tag(longbow.utils.constants.READ_BARCODE_CONF_FACTOR_TAG)}\n")

    if ignore_cbc_and_umi:
        bam_out.write(a)
        return True
    else:
        # We have to check the CBC and UMI:
        # TODO: Expand this for other annotation segments.
        cbc_ok = (not model_annotates_cbc) or bam_utils.has_cbc(a)
        umi_ok = (not model_annotates_umi) or bam_utils.has_umi(a)

        if not cbc_ok:
            logger.warning(f"Read {a.query_name} has no CBC.  Ignoring.")
            return False
        elif not umi_ok:
            logger.warning(f"Read {a.query_name} has no UMI.  Ignoring.")
            return False
        else:
            bam_out.write(a)
            return True


def create_simple_split_array_element(delim_name, end_coord, model, prev_delim_name, read, segments, start_coord, split_read_index):
    """Package an array element into an AlignedSegment from the results of simple splitting rules."""

    # Add one to end_coord because coordinates are inclusive:
    a = pysam.AlignedSegment()
    a.query_sequence = f"{read.query_sequence[start_coord:end_coord + 1]}"
    a.query_qualities = read.query_alignment_qualities[start_coord: end_coord + 1]
    a.tags = read.get_tags()
    a.flag = 4  # unmapped flag
    a.mapping_quality = 255

    # Reset read name (we need a unique name for each read that's also compatible with IsoSeq3)
    movie_name = read.query_name.split("/")[0]
    a.query_name = bam_utils.generate_read_name(movie_name, read.get_tag("zm"), split_read_index)
    a.set_tag(longbow.utils.constants.READ_ALTERED_NAME_TAG, f"{read.query_name}/{start_coord}_{end_coord}/{prev_delim_name}-{delim_name}")

    # Get our annotations for this read and modify their output coordinates so that they're relative to the length of
    # this array element / read segment:
    out_segments = []
    segments_to_annotate = []
    for s in segments:
        if start_coord <= s.start <= end_coord:
            seg_info = SegmentInfo(s.name, s.start - start_coord, s.end - start_coord)
            out_segments.append(seg_info)

            # If we have to annotate this segment, store it here for annotation later:
            if (model.annotation_segments is not None) and (s.name in model.annotation_segments.keys()):
                segments_to_annotate.append(seg_info)


    # Set our segments tag to only include the segments in this read:
    a.set_tag(
        longbow.utils.constants.SEGMENTS_TAG,
        longbow.utils.constants.SEGMENT_TAG_DELIMITER.join([s.to_tag() for s in out_segments]),
    )

    # Store tags where we'll record a refined tag instead of extracting query subsequence
    ba_tags = {}
    if a.has_tag(longbow.utils.constants.READ_INDEX_ARRAY_TAG):
        for index_tag in a.get_tag(longbow.utils.constants.READ_INDEX_ARRAY_TAG).split(","):
            tag, tag_start, tag_stop, bc = re.split("[:-]", index_tag)
            tag_start = int(tag_start) - start_coord
            tag_stop = int(tag_stop) - start_coord
            ba_tags[f'{tag}:{tag_start}-{tag_stop}'] = bc

        a.set_tag(longbow.utils.constants.READ_INDEX_ARRAY_TAG, None)

    # Annotate any segments in this array element that we have to:
    clipped_tags = set()
    for s in segments_to_annotate:
        for t in model.annotation_segments[s.name]:
            field_tag_name, pos_tag_name = t
            # First annotate the segment itself:
            ba_tag = f'{s.name}:{s.start}-{s.end}'
            seq = a.query_sequence[s.start:s.end + 1]
            if len(ba_tags) > 0:
                seq = '-' if ba_tag not in ba_tags else ba_tags[ba_tag]
            a.set_tag(field_tag_name, seq)
            clipped_tags.add(seq)

            # Next annotate the starting position:
            a.set_tag(pos_tag_name, s.start)

            # Next check if the annotation is our READ_RAW_BARCODE_TAG.
            # If so, we should annotate the confidence score as well:
            if field_tag_name == longbow.utils.constants.READ_RAW_BARCODE_TAG:
                # Get the length from the model:
                # barcode_length = list(model.adapter_dict[s.name].values())[0]
                qual_bases = a.query_qualities[s.start:s.end + 1]
                conf_factor = int(np.round(bam_utils.get_confidence_factor_raw_quals(qual_bases)))

                a.set_tag(longbow.utils.constants.READ_BARCODE_CONF_FACTOR_TAG, conf_factor)

    # Set IsoSeq3-compatible tags:
    a.set_tag(longbow.utils.constants.READ_CLIPPED_SEQS_LIST_TAG, ','.join(sorted(clipped_tags)))
    a.set_tag(longbow.utils.constants.READ_NUM_CONSENSUS_PASSES_TAG, 1)
    a.set_tag(longbow.utils.constants.READ_ZMW_NAMES_TAG, read.query_name)
    a.set_tag(longbow.utils.constants.READ_NUM_ZMWS_TAG, 1)
    a.set_tag(longbow.utils.constants.READ_TAGS_ORDER_TAG,
              f'{longbow.utils.constants.READ_RAW_BARCODE_TAG}-{longbow.utils.constants.READ_RAW_UMI_TAG}')

    # Set our tag indicating that this read is now segmented:
    a.set_tag(longbow.utils.constants.READ_IS_SEGMENTED_TAG, True)

    return a
