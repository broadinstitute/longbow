import gzip
import logging
import math
import os
import re
import sys
import time
from collections import OrderedDict, defaultdict
from functools import reduce

import click
import numpy as np
import pysam
import ssw
from matplotlib import transforms
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from polyleven import levenshtein

import longbow.utils.constants

from ..utils import bam_utils, cli_utils

logger = logging.getLogger(__name__)

PROG_NAME = "inspect"

DEFAULT_COLOR_MAP_ENTRY = "DEFAULT"


@click.command(PROG_NAME)
@click.option(
    "-r",
    "--read-names",
    type=str,
    multiple=True,
    help="read names (or file(s) of read names) to inspect",
)
@cli_utils.input_pbi
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
@cli_utils.model
@click.option(
    "--seg-score",
    is_flag=True,
    default=False,
    show_default=True,
    help="Display alignment score for annotated segments.  (--quick mode only)",
)
@click.option(
    "--max-length",
    type=int,
    default=longbow.utils.constants.DEFAULT_MAX_READ_LENGTH,
    show_default=True,
    required=False,
    help="Maximum length of a read to process.  Reads beyond this length will not be annotated.  If the input file has already been annotated, this parameter is ignored.",
)
@click.option(
    "--min-rq",
    type=float,
    default=-2,
    show_default=True,
    required=False,
    help="Minimum ccs-determined read quality for a read to be annotated.  CCS read quality range is [-1,1].  If the input file has already been annotated, this parameter is ignored.",
)
@click.option(
    "-q",
    "--quick",
    is_flag=True,
    default=False,
    show_default=True,
    help="Create quick (simplified) inspection figures.",
)
@click.option(
    "-a",
    "--annotated-bam",
    type=click.File(),
    show_default=True,
    required=False,
    help="Store annotations from a downstream BAM file so they can be displayed on reads from previous processing steps.",
)
@cli_utils.input_bam
def main(
    read_names,
    pbi,
    file_format,
    outdir,
    model,
    seg_score,
    max_length,
    min_rq,
    quick,
    annotated_bam,
    input_bam,
):
    """Inspect the classification results on specified reads."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    read_name_set = set()
    for read_name in read_names:
        if os.path.isfile(read_name):
            logger.info(f"Loading read names from {read_name}")
            with open(read_name) as fh:
                read_name_set.update(line.strip() for line in fh)
        else:
            read_name_set.add(read_name)

    # overwrite read_names with the actual list of names
    read_names = sorted(read_name_set)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    anns, name_map = None, None
    if annotated_bam is not None:
        anns, name_map = _load_annotations(annotated_bam)

    # Open our bam file:
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file:
        # Get our model:
        lb_model = bam_utils.load_model(model, bam_file)
        logger.info(f"Using {lb_model.name}: {lb_model.description}")

        logger.info(f"Figure drawing mode: {'simplified' if quick else 'extended'}")

        if seg_score and not quick:
            logger.warning(
                "Seg score selected, but not using quick mode.  Ignoring --seg-score."
            )

        # If we have read names, we should use them to inspect the file:
        if len(read_name_set) > 0:
            logger.info(f"Looking for {len(read_name_set)} read names")

            pbi = f"{input_bam}.pbi" if pbi is None else pbi
            if not os.path.exists(pbi):
                # raise FileNotFoundError(f"Missing .pbi file for {input_bam}")
                logger.info(
                    "No .pbi file available. Inspecting whole input bam file until we find specified reads."
                )
                for read in bam_file:
                    if read.query_name in read_name_set:
                        __create_read_figure(
                            file_format,
                            lb_model,
                            outdir,
                            read,
                            seg_score,
                            max_length,
                            min_rq,
                            quick,
                            anns,
                            name_map,
                        )
            else:
                file_offsets = load_read_offsets(pbi, load_read_names(read_names))

                for i, z in enumerate(file_offsets):
                    if not file_offsets[z]["offset"]:
                        logger.error("Read not in index file: %s", read_names[i])
                        sys.exit(1)

                    bam_file.seek(file_offsets[z]["offset"])
                    read = bam_file.__next__()
                    __create_read_figure(
                        file_format,
                        lb_model,
                        outdir,
                        read,
                        seg_score,
                        max_length,
                        min_rq,
                        quick,
                        anns,
                        name_map,
                    )
        else:
            # Without read names we just inspect every read in the file:
            logger.info(
                "No read names given. Inspecting every read in the input bam file."
            )
            for read in bam_file:
                __create_read_figure(
                    file_format,
                    lb_model,
                    outdir,
                    read,
                    seg_score,
                    max_length,
                    min_rq,
                    quick,
                    anns,
                    name_map,
                )

    logger.info(f"Done. Elapsed time: {time.time() - t_start:2.2f}s.")


def __create_read_figure(
    file_format,
    lb_model,
    outdir,
    read,
    seg_score,
    max_length,
    min_rq,
    quick,
    anns,
    name_map,
):
    """Create a figure for the given read."""

    out = f'{outdir}/{re.sub("/", "_", read.query_name)}.{file_format}'

    seq, path, logp = annotate_read(read, lb_model, max_length, min_rq)

    if seq is not None:
        logger.info("Drawing read '%s' to '%s'", read.query_name, out)

        if quick:
            draw_simplified_state_sequence(
                seq,
                path,
                logp,
                read,
                out,
                seg_score,
                lb_model,
                size=13,
                family="monospace",
            )
        else:
            draw_extended_state_sequence(
                seq,
                path,
                logp,
                read,
                out,
                seg_score,
                lb_model,
                anns,
                name_map,
                size=13,
                family="monospace",
            )


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
                            "Read name on line %d of file %s doesn't appear to be a PacBio read: %s",
                            i,
                            r,
                            rn,
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

    fmt = bam_utils.get_pbi_format()

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

        file_offsets_hash[f"{hole_number}_{range_end}"] = {
            "read_name": read_name,
            "offset": None,
        }

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
    fseq = read.query_sequence
    fppath = []
    flogp = -math.inf

    if read.has_tag(longbow.utils.constants.SEGMENTS_CIGAR_TAG):
        logger.info("loading annotation")
        # Set our ppath from the bam file:
        fppath = re.split(
            longbow.utils.constants.SEGMENT_TAG_DELIMITER,
            read.get_tag(longbow.utils.constants.SEGMENTS_CIGAR_TAG),
        )

        # Set our logp from the bam file:
        flogp = read.get_tag(longbow.utils.constants.READ_MODEL_SCORE_TAG)
    else:
        logger.info("annotating read")
        # Check for max length and min quality:
        if len(read.query_sequence) > max_length:
            logger.warning(
                f"Read is longer than max length.  "
                f"Skipping: {read.query_name} ({len(read.query_sequence)} > {max_length})"
            )
            return [None] * 3
        elif read.get_tag("rq") < min_rq:
            logger.warning(
                f"Read quality is below the minimum.  "
                f"Skipping: {read.query_name} ({read.get_tag('rq')} < {min_rq})"
            )
            return [None] * 3
        for seq in [
            read.query_sequence,
            bam_utils.reverse_complement(read.query_sequence),
        ]:
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
    # three_p_adapter_color = "#EFC94C"
    random_color = "#aaaaaa"

    cbc_color = "#11FF66"
    # umi_color = "#1166FF"

    # Color for states we haven't enumerated:
    default_color = "#c2a8f0"

    color_map = {
        longbow.utils.constants.RANDOM_SEGMENT_NAME: random_color,
        DEFAULT_COLOR_MAP_ENTRY: default_color,
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


def adjust_color_for_existing_colors(
    color, existing_color_strings, color_step, min_dist=50
):
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
            color_dist = np.sqrt(
                reduce(
                    lambda x, y: x + y,
                    list(map(lambda x: (x[1] - x[0]) ** 2, zip(color, other_color))),
                )
            )

            if color_dist < min_dist:
                color_dist_ok = False
                color = [(c + color_step) % 256 for c in color]

        if color_dist_ok:
            break

        loop_num += 1

    if loop_num >= max_num_loops:
        logger.warning(f"Could not get distinct enough color.  Using color: {color}")

    return color


def _adjust_aligned_state_sequence(seq, path, library_model, read_anns):
    cur_state = ""
    cur_adapter = None
    cur_pos = 0
    cur_seq = ""

    new_path = []

    last_path = ""

    seq_pos = 0
    for p in path:
        if cur_state == "CBC" or cur_state == "UMI":
            tag = "CB" if cur_state == "CBC" else "ZU"
            ssw_aligner = ssw.Aligner()
            s = None
            si = 0
            for i in range(len(read_anns)):
                s1 = ssw_aligner.align(cur_seq, read_anns[i][tag])
                if s is None:
                    s = s1
                    si = i
                else:
                    if s1.score > s.score:
                        s = s1
                        si = i

            new_cigar = []
            for opgroup in list(filter(None, re.split(r"(\d+[MIDS])", s.cigar))):
                q = re.match(r"(\d+)([MIDS])", opgroup)
                op = q.group(2)
                oplen = int(q.group(1))
                new_cigar.append(f"{op}{oplen}")

            last_path = f"{cur_state}:{''.join(new_cigar)}:{read_anns[si][tag]}"

        if last_path != "":
            new_path.append(last_path)

        state, ops = re.split(":", p)
        last_path = p

        for opgroup in list(filter(None, re.split(r"(R?[MID]A?B?\d+)", ops))):
            q = re.match(r"(R?[MID]A?B?)(\d+)", opgroup)
            op = q.group(1)
            oplen = int(q.group(2))

            if state != cur_state:
                cur_state = state
                cur_pos = 0
                cur_adapter = None
                cur_seq = ""

                if (
                    cur_state in library_model.adapter_dict
                    and cur_state not in library_model.annotation_segments
                    and type(library_model.adapter_dict[cur_state]) == str
                ):
                    cur_adapter = library_model.adapter_dict[cur_state]

            for _ in range(oplen):
                base = seq[seq_pos] if seq_pos < len(seq) else " "

                if op == "M":
                    if cur_adapter is not None:
                        cur_pos += 1

                    seq_pos += 1
                    cur_seq = f"{cur_seq}{base}"
                elif op == "D" or op == "RD":
                    if cur_adapter is not None:
                        cur_pos += 1
                elif op == "I" or op == "RI":
                    seq_pos += 1
                    cur_seq = f"{cur_seq}{base}"

    if last_path != "":
        new_path.append(last_path)

    return new_path


def _make_aligned_state_sequence(seq, path, library_model):
    observed_track = []
    mismatch_track = []
    expected_track = []
    classification_track = []

    cur_state = ""
    cur_adapter = None
    cur_pos = 0

    seq_pos = 0
    for p in path:
        f = re.split(":", p)
        state, ops = f[0], f[1]

        for opgroup in list(filter(None, re.split(r"(R?[MIDS]A?B?\d+)", ops))):
            q = re.match(r"(R?[MIDS]A?B?)(\d+)", opgroup)
            op = q.group(1)
            oplen = int(q.group(2))

            if state != cur_state:
                cur_state = state
                cur_pos = 0
                cur_adapter = None

                if (
                    cur_state in library_model.adapter_dict
                    and cur_state not in library_model.annotation_segments
                    and type(library_model.adapter_dict[cur_state]) == str
                ):
                    cur_adapter = library_model.adapter_dict[cur_state]
                elif len(f) > 2:
                    cur_adapter = f[2]

            for _ in range(oplen):
                base = seq[seq_pos] if seq_pos < len(seq) else " "

                if op == "M":
                    classification_track.append(state)
                    observed_track.append(base)

                    if cur_adapter is not None:
                        expected_track.append(cur_adapter[cur_pos])

                        if cur_adapter[cur_pos] == base:
                            mismatch_track.append("|")
                        else:
                            mismatch_track.append(".")

                        cur_pos += 1
                    else:
                        expected_track.append(" ")
                        mismatch_track.append(" ")

                    seq_pos += 1
                elif op == "D" or op == "RD":
                    classification_track.append(state)
                    observed_track.append("-" if op == "D" else " ")

                    if cur_adapter is not None:
                        expected_track.append(cur_adapter[cur_pos])
                        mismatch_track.append(" ")
                        cur_pos += 1
                    else:
                        expected_track.append(" ")
                        mismatch_track.append(" ")
                elif op == "I" or op == "RI":
                    classification_track.append(state)
                    observed_track.append(base)

                    if cur_adapter is not None:
                        expected_track.append("-" if op == "I" else " ")
                        mismatch_track.append(" ")
                    else:
                        expected_track.append(" ")
                        mismatch_track.append(" ")

                    seq_pos += 1

    return observed_track, mismatch_track, expected_track, classification_track


def format_state_sequence(seq, path, library_model, line_length=150):
    color_map = create_colormap_for_model(library_model)

    labelled_bases = []
    state_labels = []
    state_colors = []

    for i in range(0, len(path), line_length):
        bases = []
        label = re.split(":", path[i])[0]

        for j in range(i, min(i + line_length, len(path))):
            a = seq[j]
            b = re.split(":", path[j])[0]

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


def draw_extended_state_sequence(
    seq, path, logp, read, out, show_seg_score, library_model, anns, name_map, **kwargs
):
    line_length = 150

    if name_map is not None:
        path = _adjust_aligned_state_sequence(
            seq,
            path,
            library_model,
            list(map(lambda x: anns[x], name_map[read.query_name])),
        )

    (
        observed_track,
        mismatch_track,
        expected_track,
        classification_track,
    ) = _make_aligned_state_sequence(seq, path, library_model)

    f = Figure(figsize=(24, 24))
    ax = f.add_axes([0.1, 0.1, 0.8, 0.8])
    t = ax.transData

    collapsed_annotations = bam_utils.collapse_annotations(path)
    _set_figure_title(f, library_model, logp, collapsed_annotations, read)

    f.patch.set_visible(False)
    ax.axis("off")

    columns = line_length
    rows = 4 * 80

    ax.set_xlim([0, columns])
    ax.set_ylim([0, rows])

    row = 0
    column = 0
    total_bases_seen = 0

    color_map = create_colormap_for_model(library_model)

    last_label = ""
    for i, (o, m, e, c) in enumerate(
        zip(observed_track, mismatch_track, expected_track, classification_track)
    ):
        color = color_map[c]

        if column == 0:
            ax.text(
                column - 1,
                rows - row,
                total_bases_seen,
                color="#000000",
                transform=t,
                va="bottom",
                ha="right",
            )

        if i == len(observed_track) - 1:
            ax.text(
                column + 2,
                rows - row,
                len(read.query_sequence),
                color="#000000",
                transform=t,
                va="bottom",
                ha="left",
            )

        ax.text(
            column,
            rows - row,
            o,
            color=color,
            transform=t,
            bbox=dict(facecolor=f"{color}22", edgecolor=f"{color}22", pad=0),
            **kwargs,
        )

        ax.text(
            column,
            rows - row - 4,
            m,
            transform=t,
            **kwargs,
        )

        if e != " ":
            ax.text(
                column,
                rows - row - 8,
                e,
                color=color,
                transform=t,
                bbox=dict(facecolor=f"{color}22", edgecolor=f"{color}22", pad=0),
                **kwargs,
            )

        if c != last_label:
            last_label = c

            ax.text(
                column,
                rows - row + 4,
                c,
                transform=t,
                va="bottom",
                ha="left",
                fontsize=8,
                bbox=dict(facecolor="white", edgecolor="black"),
            )

        total_bases_seen += 1
        column += 1
        if column >= columns:
            column = 0
            row += 15

    FigureCanvasAgg(f).print_figure(out, bbox_inches="tight")


def _expand_cigar_sequence(cigar_path):
    """Expand a cigar sequence to a `ppath` list ala the model."""
    ppath = []
    op_re = re.compile("([A-Za-z]+)(\d+)")
    for p in cigar_path:
        # Get our segment:
        seg = re.split(":", p)[0]

        # Get our cigar operators for this segment:
        ops_string = re.split(":", p)[1]
        while len(ops_string) >= 1:
            m = op_re.match(ops_string)
            op, count = m.groups()

            # Update our ops string:
            ops_string = ops_string[m.span()[1] :]

            # We need to ignore deletions:
            if op != "D":
                # Add the correct number of operations for this cigar to the ppath:
                for i in range(int(count)):
                    ppath.append(f"{seg}:{i}")

    return ppath


def draw_simplified_state_sequence(
    seq, path, logp, read, out, show_seg_score, library_model, **kwargs
):
    line_length = 150

    # TODO: This is a bandaid.  This method should not exist.  Instead we should fix `format_state_sequence`.
    ppath = _expand_cigar_sequence(path)
    base_strings, colors, labels = format_state_sequence(
        seq, ppath, library_model, line_length=line_length
    )

    # Collapse the path here so we can later annotate our reads with the correct scores:
    segments = bam_utils.collapse_annotations(path) if show_seg_score else None

    f = Figure(figsize=(24, 24))

    ax = f.add_axes([0.1, 0.1, 0.8, 0.8])
    t = ax.transData

    collapsed_annotations = bam_utils.collapse_annotations(path)
    _set_figure_title(f, library_model, logp, collapsed_annotations, read)

    f.patch.set_visible(False)
    ax.axis("off")

    columns = line_length
    rows = 4 * 80

    ax.set_xlim([0, columns])
    ax.set_ylim([0, rows])

    row = 0
    column = 0
    row_letters_seen = 0
    total_bases_seen = 0

    last_label = None
    total_segments_seen = 0

    for base_string, color, lbl in zip(base_strings, colors, labels):
        if column == 0:
            ax.text(0, rows - row, f"{total_bases_seen}  ", transform=t, ha="right")

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
        ex = text.get_window_extent()

        # I truly have no idea why I need the scaling factor, but if I don't have this, PDF images
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
            elif (
                library_model.has_named_random_segments
                and lbl in library_model.named_random_segments
            ):
                # Handle basic named random segment:
                if (
                    library_model.adapter_dict[lbl]
                    == longbow.utils.constants.RANDOM_SEGMENT_NAME
                ):
                    length = len(collapsed_annotations[total_segments_seen])
                    seg_score_string = f" [{length}]"

                # Handle special random segments:
                elif type(library_model.adapter_dict[lbl]) is dict:
                    special_seg_type = list(library_model.adapter_dict[lbl].keys())[0]

                    if (
                        special_seg_type
                        == longbow.utils.constants.FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME
                    ):
                        # Mark the length that this known random segment should be:
                        length = list(library_model.adapter_dict[lbl].values())[0]
                        seg_score_string = f" [{length}]"
                    else:
                        logger.warning(
                            f"Ignoring score/length for unknown special segment type: {special_seg_type}"
                        )

            # If we want to show the segment scores, we calculate them here:
            elif show_seg_score and (
                not library_model.has_named_random_segments
                or lbl not in library_model.named_random_segments
            ):
                if type(library_model.adapter_dict[lbl]) is dict:
                    special_seg_type = list(library_model.adapter_dict[lbl].keys())[0]

                    if (
                        special_seg_type
                        == longbow.utils.constants.HPR_SEGMENT_TYPE_NAME
                    ):
                        base, count = library_model.adapter_dict[lbl][special_seg_type]
                        known_segment_seq = base * count
                        segment_bases = _get_segment_bases(
                            seq, total_bases_seen, segments
                        )

                        # Annotate the score
                        seg_score_string = (
                            f" ("
                            f"{len(known_segment_seq) - levenshtein(segment_bases, known_segment_seq)}"
                            f"/{len(known_segment_seq)})"
                        )

                        # Also annotate the actual length:
                        seg_score_string = f"{seg_score_string} [{len(collapsed_annotations[total_segments_seen])}]"

                    else:
                        logger.warning(
                            f"Ignoring score/length for unknown special segment type: {special_seg_type}"
                        )
                else:
                    known_segment_seq = library_model.adapter_dict[lbl]
                    segment_bases = _get_segment_bases(seq, total_bases_seen, segments)

                    seg_score_string = (
                        f" ({len(known_segment_seq) - levenshtein(segment_bases, known_segment_seq)}"
                        f"/{len(known_segment_seq)})"
                    )

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
            t = transforms.offset_copy(
                text.get_transform(), x=0, y=-30 * row, units="dots"
            )

        # Store our label so we don't write it twice:
        last_label = lbl

    FigureCanvasAgg(f).print_figure(out, bbox_inches="tight")


def _set_figure_title(f, library_model, logp, collapsed_annotations, read):
    qual_string = (
        f"[read qual: {read.get_tag('rq'):0.03f}/1.000]    "
        if read.has_tag("rq")
        else ""
    )
    np_string = f"[# Passes: {read.get_tag('np')}]    " if read.has_tag("np") else ""
    is_rc_string = (
        "[Direction: RC]"
        if read.has_tag("RC") and read.get_tag("RC") == 1
        else "[Direction: Forward]"
    )

    # Calculate whether the read has the correct array structure:
    read_mas_adapters = [
        s.name
        for s in collapsed_annotations
        if s.name in library_model.array_model["structure"]
    ]
    (
        segment_order_valid,
        key_adapters_found,
        first_key_adapter_indx,
    ) = library_model.validate_segment_order(collapsed_annotations)

    printed_array_model_name = library_model.array_model["name"].replace("_", "\\_")
    valid_library_order_string = (
        "[$\\bf{"
        + printed_array_model_name
        + "}$ adapters: "
        + f"{' '.join(library_model.key_adapters)}]"
    )
    is_valid_order_string = (
        "[Segment order: $\\bf{Valid}$]"
        if segment_order_valid
        else "[Segment order: $\\bf{INVALID}$]"
    )
    read_mas_adapter_string = f"[Key adapters: {' '.join(read_mas_adapters)}]"
    f.suptitle(
        r"$\bf{"
        + read.query_name.replace("_", "\\_")
        + "}$"
        + f"\n{qual_string}{np_string}"
        f"[{len(read.query_sequence)} bp]    {is_rc_string}    "
        r"[$\bf{"
        + library_model.name.replace("_", "\\_")
        + "}$"
        + f" model score: {logp:.2f}]\n"
        f"{is_valid_order_string}    {read_mas_adapter_string}    {valid_library_order_string}",
        fontsize=16,
    )


def _get_segment_bases(seq, position, segments):
    """Get the bases for the whole segment based on the given position in the given sequence."""
    segment_bases = ""
    for s in segments:
        if s.start <= position <= s.end:
            segment_bases = seq[s.start : s.end + 1]
            break

    return segment_bases.upper()


def _load_annotations(annotated_bam):
    anns = defaultdict(dict)
    name_map = defaultdict(set)

    pysam.set_verbosity(0)
    with pysam.AlignmentFile(
        annotated_bam, "rb", check_sq=False, require_index=False
    ) as annotated_bam_file:
        for r in annotated_bam_file:
            for k, v in r.get_tags():
                anns[r.query_name][k] = v

            name_map[r.get_tag("im")].add(r.query_name)

    return anns, name_map
