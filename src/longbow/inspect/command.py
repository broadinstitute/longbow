import logging
import sys
import os
import time
import re
import math
from functools import reduce

import click
import click_log

import gzip

import numpy as np
import pysam

import editdistance

from collections import OrderedDict
from construct import *

import matplotlib.pyplot as plt
from matplotlib import transforms

import longbow.utils.constants
from ..utils import model
from ..utils.model import LibraryModel
from ..utils import bam_utils

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("inspect")
click_log.basic_config(logger)

DEFAULT_COLOR_MAP_ENTRY = "DEFAULT"


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-r",
    "--read-names",
    type=str,
    multiple=True,
    help="read names (or file(s) of read names) to inspect",
)
@click.option(
    "-p",
    "--pbi",
    required=False,
    type=click.Path(exists=True),
    help="BAM .pbi index file",
)
@click.option(
    "-f",
    "--file-format",
    default="pdf",
    show_default=True,
    type=click.Choice(["png", "pdf"]),
    help="Image file format",
)
@click.option(
    "-o",
    "--outdir",
    default=".",
    show_default=True,
    required=False,
    type=click.Path(exists=False),
    help="Output directory",
)
@click.option(
    "-m",
    "--model",
    help="The model to use for annotation.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
)
@click.option(
    '--seg-score',
    is_flag=True,
    default=False,
    show_default=True,
    help="Display alignment score for annotated segments."
)
@click.option(
    "--max-length",
    type=int,
    default=longbow.utils.constants.DEFAULT_MAX_READ_LENGTH,
    show_default=True,
    required=False,
    help="Maximum length of a read to process.  Reads beyond this length will not be annotated.  If the input file has already been annotated, this parameter is ignored."
)
@click.option(
    "--min-rq",
    type=float,
    default=-2,
    show_default=True,
    required=False,
    help="Minimum ccs-determined read quality for a read to be annotated.  CCS read quality range is [-1,1].  If the input file has already been annotated, this parameter is ignored."
)
@click.argument("input-bam", type=click.Path(exists=True))
def main(read_names, pbi, file_format, outdir, model, seg_score, max_length, min_rq, input_bam):
    """Inspect the classification results on specified reads."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Open our bam file:
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file:
        # Get our model:
        if model is None:
            lb_model = LibraryModel.from_json_obj(bam_utils.get_model_from_bam_header(bam_file.header))
        elif model is not None and LibraryModel.has_prebuilt_model(model):
            lb_model = LibraryModel.build_pre_configured_model(model)
        else:
            lb_model = LibraryModel.from_json_file(model)

        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        # If we have read names, we should use them to inspect the file:
        if len(read_names) > 0:

            pbi = f"{input_bam}.pbi" if pbi is None else pbi
            if not os.path.exists(pbi):
                # raise FileNotFoundError(f"Missing .pbi file for {input_bam}")
                logger.info("No .pbi file available. Inspecting whole input bam file until we find specified reads.")
                for read in bam_file:
                    if read.query_name in read_names:
                        __create_read_figure(file_format, lb_model, outdir, read, seg_score, max_length, min_rq)
            else:
                file_offsets = load_read_offsets(pbi, load_read_names(read_names))

                for i, z in enumerate(file_offsets):

                    if not file_offsets[z]["offset"]:
                        logger.error("Read not in index file: %s", read_names[i])
                        sys.exit(1)

                    bam_file.seek(file_offsets[z]["offset"])
                    read = bam_file.__next__()
                    __create_read_figure(file_format, lb_model, outdir, read, seg_score, max_length, min_rq)
        else:
            # Without read names we just inspect every read in the file:
            logger.info("No read names given. Inspecting every read in the input bam file.")
            for read in bam_file:
                __create_read_figure(file_format, lb_model, outdir, read, seg_score, max_length, min_rq)

    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)


def __create_read_figure(file_format, lb_model, outdir, read, seg_score, max_length, min_rq):
    """Create a figure for the given read."""

    out = f'{outdir}/{re.sub("/", "_", read.query_name)}.{file_format}'

    seq, path, logp = annotate_read(read, lb_model, max_length, min_rq)

    if seq is not None:
        logger.info("Drawing read '%s' to '%s'", read.query_name, out)
        draw_state_sequence(seq, path, logp, read, out, seg_score, lb_model, size=13, family="monospace")


def load_read_names(read_name_args):
    read_names = []

    for r in read_name_args:
        if os.path.exists(r):
            with open(r, "r") as f:
                for i, line in enumerate(f):
                    # Space / Tab / Newline / Line feed are all forbidden in read names by the sam spec, so we can
                    # trim it all off:
                    rn = line.strip()
                    if not bam_utils.PB_READ_NAME_RE.match(rn):
                        logger.error(
                            "Read name on line %d of file %s doesn't appear to be a PacBio read: %s", i, r, rn
                        )
                        sys.exit(1)

                    read_names.append(rn)
        else:
            if not bam_utils.PB_READ_NAME_RE.match(r):
                logger.error("Read name doesn't appear to be a PacBio read: %s", r)
                sys.exit(1)
            read_names.append(r)

    return read_names


def load_read_offsets(pbi_file, read_names):
    """
    Compute file offsets for specified read names
    """

    # Decode PacBio .pbi file.  This is not a full decode of the index, only the parts we need for read selection.
    # More on index format at https://pacbiofileformats.readthedocs.io/en/9.0/PacBioBamIndex.html .

    fmt = Struct(
        # Header
        "magic" / Const(b"PBI\x01"),
        "version_patch" / Int8ul,
        "version_minor" / Int8ul,
        "version_major" / Int8ul,
        "version_empty" / Int8ul,
        "pbi_flags" / Int16ul,
        "n_reads" / Int32ul,
        "reserved" / Padding(18),
        # Basic information section (columnar format)
        "rgId" / Padding(this.n_reads * 4),
        "qStart" / Padding(this.n_reads * 4),
        "qEnd" / Array(this.n_reads, Int32sl),
        "holeNumber" / Array(this.n_reads, Int32sl),
        "readQual" / Padding(this.n_reads * 4),
        "ctxtFlag" / Padding(this.n_reads * 1),
        "fileOffset" / Array(this.n_reads, Int64sl),
    )

    # Make a list of bgzf virtual file offsets for sharding and store ZMW counts.
    file_offsets_hash = OrderedDict()
    for read_name in read_names:
        name_pieces = re.split("/", read_name)

        hole_number = int(name_pieces[1])
        if "ccs" in name_pieces:
            if len(name_pieces) > 3:
                rr = re.split("_", name_pieces[3])
                range_end = int(rr[1]) - int(rr[0]) + 1
            else:
                range_end = -1
        else:
            range_end = re.split("_", name_pieces[2])[1]

        file_offsets_hash[f"{hole_number}_{range_end}"] = {"read_name": read_name, "offset": None}

    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        for j in range(0, idx_contents.n_reads):
            key1 = f"{idx_contents.holeNumber[j]}_{idx_contents.qEnd[j]}"
            key2 = f"{idx_contents.holeNumber[j]}_-1"

            # Save virtual file offsets for the selected reads only
            if key1 in file_offsets_hash:
                file_offsets_hash[key1]["offset"] = idx_contents.fileOffset[j]
            elif key2 in file_offsets_hash:
                file_offsets_hash[key2]["offset"] = idx_contents.fileOffset[j]

    return file_offsets_hash


def annotate_read(read, m, max_length, min_rq):
    flogp = -math.inf
    fseq = read.query_sequence
    fppath = []

    if read.has_tag(longbow.utils.constants.SEGMENTS_TAG):
        tag = re.split(",", read.get_tag(longbow.utils.constants.SEGMENTS_TAG))

        for e in tag:
            state, rrange = re.split(":", e)
            qStart, qEnd = re.split("-", rrange)
            qLen = int(qEnd) - int(qStart) + 1

            fppath.extend([state] * qLen)

        # Set our logp from the bam file:
        flogp = read.get_tag(longbow.utils.constants.READ_MODEL_SCORE_TAG)
    else:

        # Check for max length and min quality:
        if len(read.query_sequence) > max_length:
            logger.warning(f"Read is longer than max length.  "
                           f"Skipping: {read.query_name} ({len(read.query_sequence)} > {max_length})")
            return [None] * 3
        elif read.get_tag("rq") < min_rq:
            logger.warning(f"Read quality is below the minimum.  "
                           f"Skipping: {read.query_name} ({read.get_tag('rq')} < {min_rq})")
            return [None] * 3
        for seq in [read.query_sequence, bam_utils.reverse_complement(read.query_sequence)]:
            logp, ppath = m.annotate(seq)

            if logp > flogp:
                fseq = seq
                flogp = logp
                fppath = ppath

    return fseq, fppath, flogp


def create_colormap_for_model(m):
    adapter_state_color = "#DF5A49"
    scaffold_state_color = "#334D5C"
    poly_a_color = "#45B29D"
    coding_region_color = "#4455FF"
    three_p_adapter_color = "#EFC94C"
    random_color = "#aaaaaa"

    cbc_color = "#11FF66"
    umi_color = "#1166FF"

    # Color for states we haven't enumerated:
    default_color = "#c2a8f0"

    color_map = {
        longbow.utils.constants.RANDOM_SEGMENT_NAME: random_color,
        DEFAULT_COLOR_MAP_ENTRY: default_color
    }

    # Generate the rest of the colors:
    for name in m.get_all_node_names():
        if name in color_map:
            continue
        if len(name) == 1:
            c = adapter_state_color
        elif m.has_coding_region and name == m.coding_region:
            c = coding_region_color
        elif name in longbow.utils.constants.MAS_SCAFFOLD_NAMES:
            c = scaffold_state_color
        elif m.has_named_random_segments and name in m.named_random_segments:
            # We'll do these next:
            continue
        elif name.upper().startswith("POLY"):
            c = poly_a_color
        else:
            c = default_color

        color_map[name] = c

    if m.has_named_random_segments:
        # Now generate random segment colors:
        # Start from the color of the coding region:
        rc = [
            int(cbc_color[1:3], 16),
            int(cbc_color[3:5], 16),
            int(cbc_color[5:7], 16),
        ]
        step = int(256 / len(m.named_random_segments))
        for name in m.named_random_segments:
            if m.has_coding_region and name == m.coding_region:
                continue
            c = [(c + step) % 256 for c in rc]
            c = adjust_color_for_existing_colors(c, color_map.values(), step)

            color_map[name] = f"#{c[0]:02X}{c[1]:02X}{c[2]:02X}"
            rc = c

    return color_map


def adjust_color_for_existing_colors(color, existing_color_strings, color_step, min_dist=50):
    max_num_loops = 10
    loop_num = 0

    while loop_num < max_num_loops:
        color_dist_ok = True
        for color_string in existing_color_strings:
            other_color = [
                int(color_string[1:3], 16),
                int(color_string[3:5], 16),
                int(color_string[5:7], 16),
            ]

            # Simple Euclidian distance in color space:
            color_dist = np.sqrt(reduce(lambda x, y: x+y, list(map(lambda x: (x[1]-x[0])**2, zip(color, other_color)))))

            if color_dist < min_dist:
                color_dist_ok = False
                color = [(c + color_step) % 256 for c in color]

        if color_dist_ok:
            break

        loop_num += 1

    if loop_num >= max_num_loops:
        logger.warning(f"Could not get distinct enough color.  Using color: {color}")

    return color


def format_state_sequence(seq, path, library_model, line_length=150):
    color_map = create_colormap_for_model(library_model)

    labelled_bases = []
    state_labels = []
    state_colors = []

    for i in range(0, len(path), line_length):
        bases = []
        label = path[i]

        for j in range(i, min(i + line_length, len(path))):
            a = seq[j]
            b = path[j]

            # Handle the line breaks:
            if label == b:
                bases.append(a)
            elif label != b:
                try:
                    color = color_map[label]
                except KeyError:
                    color = color_map[DEFAULT_COLOR_MAP_ENTRY]

                labelled_bases.append("".join(bases))
                state_labels.append(label)
                state_colors.append(color)

                bases = [a]
                label = b

        labelled_bases.append("".join(bases))
        state_labels.append(label)
        try:
            color = color_map[label]
        except KeyError:
            color = color_map[DEFAULT_COLOR_MAP_ENTRY]
        state_colors.append(color)

    return labelled_bases, state_colors, state_labels


def draw_state_sequence(seq, path, logp, read, out, show_seg_score, library_model, **kwargs):

    line_length = 150

    base_strings, colors, labels = format_state_sequence(seq, path, library_model, line_length=line_length)

    # Collapse the path here so we can later annotate our reads with the correct scores:
    segments = bam_utils.collapse_annotations(path) if show_seg_score else None

    f = plt.figure(figsize=(24, 24))

    ax = plt.gca()
    t = ax.transData
    canvas = ax.figure.canvas

    qual_string = f"[read qual: {read.get_tag('rq'):0.03f}/1.000]    " if read.has_tag("rq") else ""
    np_string = f"[# Passes: {read.get_tag('np')}]    " if read.has_tag("np") else ""
    is_rc_string = "[Direction: RC]" if read.has_tag("RC") and read.get_tag("RC") == 1 else "[Direction: Forward]"

    # Calculate whether the read has the correct array structure:
    collapsed_annotations = bam_utils.collapse_annotations(path)
    read_mas_adapters = [s.name for s in collapsed_annotations if len(s.name) == 1]
    segment_order_valid, key_adapters_found, first_key_adapter_indx = \
        library_model.validate_segment_order([s.name for s in collapsed_annotations])

    valid_library_order_string = f"[{library_model.name} adapters: {' '.join(library_model.key_adapters)}]"
    is_valid_order_string = "[Segment order: Valid]" if segment_order_valid else "[Segment order: INVALID]"
    read_mas_adapter_string = f"[Key adapters: {' '.join(read_mas_adapters)}]"

    f.suptitle(
        r"$\bf{" + read.query_name.replace("_", "\\_") + "}$" + f"\n{qual_string}{np_string}"
        f"[{len(read.query_sequence)} bp]    {is_rc_string}    "
        r"[$\bf{" + library_model.name.replace("_", "\\_") + "}$" + f" model score: {logp:.2f}]\n"
        f"{is_valid_order_string}    {read_mas_adapter_string}    {valid_library_order_string}",
        fontsize=16
    )

    f.patch.set_visible(False)
    ax.axis("off")

    columns = line_length
    rows = 4 * 80

    plt.xlim([0, columns])
    plt.ylim([0, rows])

    row = 0
    column = 0
    row_letters_seen = 0
    total_bases_seen = 0

    last_label = None
    total_segments_seen = 0

    for base_string, color, lbl in zip(base_strings, colors, labels):
        if column == 0:
            ytic = ax.text(0, rows - row, f"{total_bases_seen}  ", transform=t, ha="right")

        # Write classified sequence
        text = ax.text(
            0,
            rows - row,
            base_string,
            color=color,
            transform=t,
            bbox=dict(facecolor=f"{color}22", edgecolor=f"{color}22", pad=0),
            **kwargs,
        )
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()

        # I truly have no idea why I need the 0.72 scaling factor, but if I don't have this, PDF images
        # are super broken.
        scaling_factor = 1.00 if out.endswith(".png") else 0.72
        t = transforms.offset_copy(
            text.get_transform(), x=scaling_factor * ex.width, units="dots"
        )

        # Write state label if we haven't already written it:
        if lbl != last_label:
            seg_score_string = ""
            if lbl == longbow.utils.constants.RANDOM_SEGMENT_NAME:
                # Let's add the length of the random segment to the label:
                length = len(collapsed_annotations[total_segments_seen])
                seg_score_string = f" [{length}]"

            # Always display the length of known random segments of fixed length:
            elif library_model.has_named_random_segments and lbl in library_model.named_random_segments:

                # Handle basic named random segment:
                if library_model.adapter_dict[lbl] == longbow.utils.constants.RANDOM_SEGMENT_NAME:
                    length = len(collapsed_annotations[total_segments_seen])
                    seg_score_string = f" [{length}]"

                # Handle special random segments:
                elif type(library_model.adapter_dict[lbl]) is dict:

                    special_seg_type = list(library_model.adapter_dict[lbl].keys())[0]

                    if special_seg_type == longbow.utils.constants.FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                        # Mark the length that this known random segment should be:
                        length = list(library_model.adapter_dict[lbl].values())[0]
                        seg_score_string = f" [{length}]"
                    else:
                        logger.warning(f"Ignoring score/length for unknown special segment type: {special_seg_type}")

            # If we want to show the segment scores, we calculate them here:
            elif show_seg_score and (not library_model.has_named_random_segments or lbl not in library_model.named_random_segments):

                if type(library_model.adapter_dict[lbl]) is dict:
                    special_seg_type = list(library_model.adapter_dict[lbl].keys())[0]

                    if special_seg_type == longbow.utils.constants.HPR_SEGMENT_TYPE_NAME:
                        base, count = library_model.adapter_dict[lbl][special_seg_type]
                        known_segment_seq = base * count
                        segment_bases = _get_segment_bases(seq, total_bases_seen, segments)

                        # Annotate the score
                        seg_score_string = f" ({len(known_segment_seq) - editdistance.eval(segment_bases, known_segment_seq)}" \
                                           f"/{len(known_segment_seq)})"

                        # Also annotate the actual length:
                        seg_score_string = f"{seg_score_string} [{len(collapsed_annotations[total_segments_seen])}]"

                    else:
                        logger.warning(f"Ignoring score/length for unknown special segment type: {special_seg_type}")
                else:
                    known_segment_seq = library_model.adapter_dict[lbl]
                    segment_bases = _get_segment_bases(seq, total_bases_seen, segments)

                    seg_score_string = f" ({len(known_segment_seq) - editdistance.eval(segment_bases, known_segment_seq)}" \
                                       f"/{len(known_segment_seq)})"

            ax.text(
                0 - (len(base_string) / 2),
                rows - row - 3.5,
                f"{lbl}{seg_score_string}",
                transform=t,
                va="bottom",
                ha="center",
                fontsize=8,
                bbox=dict(facecolor="white", edgecolor="black"),
            )

            total_segments_seen += 1

        # TODO: Replace all references to "5p_Adapter" with references to the models themselves.
        # Add in hashes after the leading structural adapter so we can visually more easily inspect the results:
        if lbl == "5p_Adapter":
            for o in range(10, 50, 10):
                if row_letters_seen + len(base_string) + o < line_length:
                    ax.text(0, rows - row + 1.7, f"{' ' * o}'", transform=t, **kwargs)

        # Decide whether we need to break into a new row
        row_letters_seen += len(base_string)
        total_bases_seen += len(base_string)
        column += 1

        if row_letters_seen >= columns:
            row_letters_seen = 0

            row += 1
            column = 0

            text = ax.text(column, row, "")
            text.draw(canvas.get_renderer())
            t = transforms.offset_copy(
                text.get_transform(), x=0, y=-30 * row, units="dots"
            )

        # Store our label so we don't write it twice:
        last_label = lbl

    y2tic = ax.text(0, rows - row, f"  {total_bases_seen}", transform=t, ha="left")

    plt.savefig(out, bbox_inches="tight")

    plt.close()


def _get_segment_bases(seq, position, segments):
    """Get the bases for the whole segment based on the given position in the given sequence."""
    segment_bases = ""
    for s in segments:
        if s.start <= position <= s.end:
            segment_bases = seq[s.start:s.end + 1]
            break

    return segment_bases.upper()

