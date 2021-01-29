import logging
import click
import click_log

import gzip
from collections import OrderedDict
from construct import *

import matplotlib.pyplot as plt
from matplotlib import transforms

from ..utils.model import *

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command(name="inspect")
@click_log.simple_verbosity_option(logger)
@click.option("-r", "--read-names", type=str, multiple=True, help="read names (or file(s) of read names) to inspect")
@click.option("-m", "--model", required=False, type=click.Path(exists=True), help="pre-trained model to apply")
@click.option("-p", "--pbi", required=False, type=click.Path(exists=True), help="BAM .pbi index file")
@click.option("-o", "--outdir", default=".", required=False, type=click.Path(exists=False), help="Output directory")
@click.argument('input-bam', type=click.Path(exists=True))
def main(read_names, model, pbi, outdir, input_bam):
    """Inspect the classification results on specified reads"""
    logger.info(f"annmas: inspect started")

    pbi = f'{input_bam}.pbi' if pbi is None else pbi
    if not os.path.exists(pbi):
        raise FileNotFoundError(f"Missing .pbi file for {input_bam}")

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    m = build_default_model()
    if model is not None:
        m.from_yaml(model)
        logger.info(f"Using pretrained annotation model {model}")
    else:
        logger.info(f"Using default annotation model")

    file_offsets = load_read_offsets(pbi, load_read_names(read_names))

    pysam.set_verbosity(0)
    bf = pysam.Samfile(input_bam, 'rb', check_sq=False, require_index=False)

    for z in file_offsets:
        bf.seek(file_offsets[z]['offset'])
        read = bf.__next__()

        out = f'{outdir}/{re.sub("/", "_", read.query_name)}.png'
        seq, path, logp = annotate_read(read, m)

        logger.info("Drawing read '%s' to '%s'", read.query_name, out)
        draw_state_sequence(seq, path, read.query_name, out)

    bf.close()

    logger.info("annmas: inspect finished")


def load_read_names(read_names):
    rn = []

    for r in read_names:
        if os.path.exists(r):
            [rn.append(line.rstrip('\n')) for line in open(r)]
        else:
            rn.append(r)

    return rn


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
        "qStart" / Array(this.n_reads, Int32sl),
        "qEnd" / Padding(this.n_reads * 4),
        "holeNumber" / Array(this.n_reads, Int32sl),
        "readQual" / Padding(this.n_reads * 4),
        "ctxtFlag" / Padding(this.n_reads * 1),
        "fileOffset" / Array(this.n_reads, Int64sl),
        )

    # Make a list of bgzf virtual file offsets for sharding and store ZMW counts.
    file_offsets_hash = OrderedDict()
    for read_name in read_names:
        movie_name, zmw, r = re.split("/", read_name)
        if r == 'ccs':
            r = 0
        else:
            r = re.split("_", r)[0]

        file_offsets_hash[f'{zmw}_{r}'] = {'read_name': read_name, 'offset': None}

    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        for j in range(0, idx_contents.n_reads):
            key = f'{idx_contents.holeNumber[j]}_{idx_contents.qStart[j]}'

            # Save virtual file offsets for the selected reads only
            if key in file_offsets_hash:
                file_offsets_hash[key]['offset'] = idx_contents.fileOffset[j]

    return file_offsets_hash


def annotate_read(read, m):
    flogp = -math.inf
    fseq = None
    for seq in [read.query_sequence, reverse_complement(read.query_sequence)]:
        logp, ppath = annotate(m, seq)

        if logp > flogp:
            fseq = seq
            flogp = logp
            fppath = ppath

    return fseq, fppath, flogp


def format_state_sequence(seq, path):
    color_hash = {
        "10x_Adapter": "blue",
        "5p_TSO": "blue",
        "Poly_A": "green",
        "3p_Adapter": "gold",
        "A": "red",
        "B": "red",
        "C": "red",
        "D": "red",
        "E": "red",
        "F": "red",
        "G": "red",
        "H": "red",
        "I": "red",
        "J": "red",
        "K": "red",
        "L": "red",
        "M": "red",
        "N": "red",
        "O": "red",
        "P": "red",
        "random": "darkgrey"
    }

    labelled_bases = []
    state_labels = []
    state_colors = []

    bases = []
    label = path[0]
    al = 0
    for i, (a, b) in enumerate(zip(list(seq), path)):
        if label == b:
            if al + len(a) < 150:
                bases.append(a)
                label = b
                al += len(a)
            else:
                labelled_bases.append("".join(bases))
                state_labels.append(label)
                state_colors.append(color_hash[label])

                bases = [a]
                label = b
                al = 0
        else:
            labelled_bases.append("".join(bases))
            state_labels.append(label)
            state_colors.append(color_hash[label])

            bases = [a]
            label = b

    return labelled_bases, state_colors, state_labels


def draw_state_sequence(seq, path, read_name, out, **kwargs):
    strings, colors, labels = format_state_sequence(seq, path)

    f = plt.figure(figsize=(24, 24))

    ax = plt.gca()
    t = ax.transData
    canvas = ax.figure.canvas

    f.suptitle(read_name, fontsize=16)
    f.patch.set_visible(False)
    ax.axis('off')

    columns = 150
    rows = 4 * 80

    plt.xlim([0, columns])
    plt.ylim([0, rows])

    letters_seen = 0
    row = 0
    column = 0
    for s, c, l in zip(strings, colors, labels):
        text = ax.text(0.5 + column * 1.2, rows - row, s, color=c, transform=t,
                       bbox=dict(facecolor='none', edgecolor=c), **kwargs)

        # Write classified sequence
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()
        t = transforms.offset_copy(
            text.get_transform(), x=ex.width, units='dots')

        # Write class label
        ax.text(0.5 + column * 1.3 - (len(s) / 2), rows - row - 3, l, transform=t, va='bottom', ha='center', fontsize=8,
                bbox=dict(facecolor='white', edgecolor='black'))

        # Decide whether we need to break into a new row
        letters_seen += len(s)
        column += 1

        if letters_seen >= columns:
            letters_seen = 0

            row += 1
            column = 0

            text = ax.text(column, row, "")
            text.draw(canvas.get_renderer())
            t = transforms.offset_copy(
                text.get_transform(), x=0, y=-30 * row, units='dots')

    plt.savefig(out, bbox_inches='tight')

    plt.close()

