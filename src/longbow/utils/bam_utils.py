import json
import math
import gzip
import logging
import click_log
import re
from numpy import isin
import pysam
import sys
import array
import collections
import os
import operator

from functools import reduce
from collections import OrderedDict
from construct import *
from inspect import currentframe, getdoc

from ..meta import VERSION

from .constants import RANDOM_SEGMENT_NAME, HPR_SEGMENT_TYPE_NAME, SEGMENTS_TAG, SEGMENTS_QUAL_TAG, SEGMENTS_RC_TAG, \
    SEGMENT_TAG_DELIMITER, READ_MODEL_NAME_TAG, READ_MODEL_SCORE_TAG, READ_APPROX_QUAL_TAG, READ_RAW_UMI_TAG, \
    READ_RAW_BARCODE_TAG, CONF_FACTOR_SCALE

from ..utils.model import LibraryModel
from ..utils import bam_utils


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("bam_utils")
click_log.basic_config(logger)

PB_READ_NAME_RE = re.compile("m[0-9]+e?_[0-9]{6}_[0-9]{6}/[0-9]+/.*")


# Named tuple to store alignment information:
class SegmentInfo(collections.namedtuple("SegmentInfo", ["name", "start", "end"])):

    _tag_regex = re.compile(r"(.*?):(\d+)-(\d+)")

    def __len__(self):
        # Must add 1 because the positions are inclusive coordinates:
        return self.end - self.start + 1

    def __str__(self):
        return f"SegmentInfo({self.to_tag()})"

    def to_tag(self):
        return f"{self.name}:{self.start}-{self.end}"

    @classmethod
    def from_tag(cls, tag_string):
        match = cls._tag_regex.match(tag_string)
        return SegmentInfo(match[1], int(match[2]), int(match[3]))


def load_read_count(pbi_file):
    """Compute file offsets for specified read names"""

    # Decode PacBio .pbi file.  This is not a full decode of the index, only the parts we need
    # until we get to the read count.
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
        )

    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        return idx_contents.n_reads


def compute_shard_offsets(pbi_file, num_shards):
    """
    Compute all possible shard offsets (keeping adjacent reads from the ZMW together)
    """

    # Decode PacBio .pbi file.  This is not a full decode of the index, only the parts we need for sharding.
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
        "qEnd" / Padding(this.n_reads * 4),
        "holeNumber" / Array(this.n_reads, Int32sl),
        "readQual" / Padding(this.n_reads * 4),
        "ctxtFlag" / Padding(this.n_reads * 1),
        "fileOffset" / Array(this.n_reads, Int64sl),
        )

    # Make a list of bgzf virtual file offsets for sharding and store ZMW counts.
    file_offsets_hash = OrderedDict()
    last_offset = 0
    zmw_count_hash = dict()
    with gzip.open(pbi_file, "rb") as f:
        idx_contents = fmt.parse_stream(f)

        for j in range(0, idx_contents.n_reads):
            # Save only the virtual file offset for the first ZMW hole number, so
            # that shard boundaries always keep reads from the same ZMW together.
            if idx_contents.holeNumber[j] not in file_offsets_hash:
                file_offsets_hash[idx_contents.holeNumber[j]] = idx_contents.fileOffset[j]

            last_offset = idx_contents.fileOffset[j]

            try:
                zmw_count_hash[idx_contents.holeNumber[j]] += 1
            except KeyError:
                zmw_count_hash[idx_contents.holeNumber[j]] = 1

    file_offsets = list(file_offsets_hash.values())
    shard_offsets = []
    read_counts = []
    read_nums = []

    read_num = 1
    by = int(math.ceil(len(file_offsets) / num_shards))
    for i in range(0, len(file_offsets), by):
        shard_offsets.append(file_offsets[i])
        read_counts.append(len(file_offsets[i:i + by]))
        read_nums.append(read_num)
        read_num += read_counts[-1]

    # For the last read in the file, pad the offset so the final comparison in write_shard() retains the final read.
    shard_offsets.append(file_offsets[-1] + 1)

    return shard_offsets, zmw_count_hash, idx_contents.n_reads, read_counts, read_nums


def create_bam_header_with_program_group(command_name, base_bam_header, description=None, models=None):
    """Create a pysam.AlignmentHeader object with program group (PG) information populated by the given arguments.

    This function is intended to be called from the 'main' function of a longbow subcommand because it uses reflection
    to pull in the first line of the docstring from the main function as the description (DS field)."""

    bam_header_dict = base_bam_header.to_dict()

    if not description:
        prev_frame = currentframe().f_back
        description = getdoc(prev_frame.f_globals['main']).split("\n")[0]

    # If we have a model here, we should add the description of the model to our program group:
    if models:
        description = description + "  MODEL(s): " + ", ".join([m.to_json(indent=None) for m in models])

    # Add our program group to it:
    pg_dict = {
        "ID": f"longbow-{command_name}-{VERSION}",
        "PN": "longbow",
        "VN": f"{VERSION}",
        # Use reflection to get the first line of the doc string the caller - the main function for our header:
        "DS": description,
        "CL": " ".join(sys.argv),
    }
    if "PG" in bam_header_dict:
        bam_header_dict["PG"].append(pg_dict)
    else:
        bam_header_dict["PG"] = [pg_dict]
    out_header = pysam.AlignmentHeader.from_dict(bam_header_dict)

    return out_header


def check_for_preexisting_files(file_list, exist_ok=False):
    """Checks if the files in the given file_list exist.
    If any file exists and exist_ok is False, this will exit the program.
    If any file exists and exist_ok is True, the program will continue.
    """

    # Allow users to be a little lazy with what input types they give:
    if not isinstance(file_list, list) and not isinstance(file_list, set):
        file_list = [file_list]

    do_files_exist = False
    for f in file_list:
        if os.path.exists(f) and not f == "/dev/null":
            if exist_ok:
                logger.warning(f"Output file exists: {f}.  Overwriting.")
            else:
                logger.error(f"Output file already exists: {f}!")
                do_files_exist = True
    if do_files_exist:
        sys.exit(1)


def load_models(models, input_bam=None):
    """Load LibraryModel objects from a BAM file header, built-in models, or an external JSON file."""

    lb_models = []
    sources = []
    if models is None and input_bam is not None:
        pysam.set_verbosity(0)
        with pysam.AlignmentFile(input_bam.name, "rb", check_sq=False, require_index=False) as bam_file:
            if bam_utils.bam_header_has_model(bam_file.header):
                for model_json in get_models_from_bam_header(bam_file.header):
                    lb_models.append(LibraryModel.from_json_obj(model_json))
                    sources.append("BAM header")
            else:
                logger.fatal("Model not specified and no model present in BAM header.")
                sys.exit(1)
    else:
        for model in models:
            if LibraryModel.has_prebuilt_model(model):
                lb_models.append(LibraryModel.build_pre_configured_model(model))
                sources.append("pre-built models")
            else:
                lb_models.append(LibraryModel.from_json_file(model))
                sources.append("user-supplied JSON")

    for lb_model, source in zip(lb_models, sources):
        logger.info(f"Using {lb_model.name} ({lb_model.description}) from {source}")

    return lb_models


def get_segment_score(read_sequence, segment, library_model, ssw_aligner=None):
    """Get the alignment score of the given segment against the read sequence."""

    # TODO: FIX THIS METHOD WITH NEW SCORING MODEL FROM INSPECT!
    # It should be something like:
    #                         base, count = library_model.adapter_dict[lbl][special_seg_type]
    #                         known_segment_seq = base * count
    #                         segment_bases = _get_segment_bases(seq, total_bases_seen, segments)
    #
    #                         seg_score_string = f" ({len(known_segment_seq) - editdistance.eval(segment_bases, known_segment_seq)}" \
    #                                            f"/{len(known_segment_seq)})"


    # We don't score random segments:
    if segment.name == RANDOM_SEGMENT_NAME or \
            (library_model.has_named_random_segments and segment.name in library_model.named_random_segments):
        return 0, 0

    # Create a default aligner if we weren't given one:
    if not ssw_aligner:
        ssw_aligner = ssw.Aligner()

    # get our model segment so we can reference it later:
    if type(library_model.adapter_dict[segment.name]) == str:
        model_seg_sequence = library_model.adapter_dict[segment.name]
    elif type(library_model.adapter_dict[segment.name]) == dict:
        if list(library_model.adapter_dict[segment.name].keys())[0] == HPR_SEGMENT_TYPE_NAME:
            b, l = list(library_model.adapter_dict[segment.name].values())[0]
            model_seg_sequence = b * l
        else:
            raise RuntimeError(f"Cannot calculate score for dict segment type: "
                               f"{type(library_model.adapter_dict[segment.name])}: {segment}")
    else:
        raise RuntimeError(
            f"Unknown segment type: {type(library_model.adapter_dict[segment.name])} for segment: {segment}"
        )

    # Get our alignment and our score:
    if segment.end - segment.start > 1:
        alignment = ssw_aligner.align(read_sequence[segment.start:segment.end], model_seg_sequence)
        optimal_score = alignment.score
    else:
        optimal_score = 0

    # The max score is the match score * the length of the reference segment
    max_score = len(model_seg_sequence) * ssw_aligner.matrix.get_match()

    return optimal_score, max_score


def collapse_annotations(path):
    """Collapses given path into a list of SegmentInfo objects."""
    last = ""
    start = 0
    segments = []
    i = 0
    for i, seg in enumerate(path):
        if seg != last:
            if i != 0:
                segments.append(SegmentInfo(last, start, i - 1))
            last = seg
            start = i
    # Don't forget the last one:
    segments.append(SegmentInfo(last, start, i))

    return segments


def write_annotated_read(read, segments, is_rc, logp, model, ssw_aligner, out_bam_file):
    """Write the given pysam.AlignedSegment read object to the given file with the given metadata."""

    # Obligatory log message:
    logger.debug(
        "Path for read %s (%2.2f)%s: %s",
        read.query_name,
        logp,
        " (RC)" if is_rc else "",
        segments,
    )

    # Set our tag and write out the read to the annotated file:
    read.set_tag(SEGMENTS_TAG, SEGMENT_TAG_DELIMITER.join([s.to_tag() for s in segments]))

    # Set the model info tags:
    read.set_tag(READ_MODEL_SCORE_TAG, logp)
    read.set_tag(READ_MODEL_NAME_TAG, model.name)

    # If we're reverse complemented, we make it easy and just reverse complement the read and add a tag saying
    # that the read was RC:
    read.set_tag(SEGMENTS_RC_TAG, is_rc)
    if is_rc:
        quals = read.query_qualities[::-1]
        seq = reverse_complement(read.query_sequence)
        read.query_sequence = seq
        read.query_qualities = quals

    # Get our segment scores and set them:
    total_score = 0
    total_max_score = 0
    score_strings = []
    for s in segments:
        score, max_score = get_segment_score(read.query_sequence, s, model, ssw_aligner)
        score_strings.append(f"{score}/{max_score}")
        total_score += score
        total_max_score += max_score

    read.set_tag(SEGMENTS_QUAL_TAG, SEGMENT_TAG_DELIMITER.join(score_strings))
    if total_max_score != 0:
        read.set_tag(READ_APPROX_QUAL_TAG, f"{total_score / total_max_score:.4f}")
    else:
        read.set_tag(READ_APPROX_QUAL_TAG, f"0.0")

    out_bam_file.write(read)


# IUPAC RC's from: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
# and https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
RC_BASE_MAP = {
    "N": "N",
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "Y": "R",
    "R": "Y",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "n": "n",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "y": "r",
    "r": "y",
    "s": "s",
    "w": "w",
    "k": "m",
    "m": "k",
    "b": "v",
    "v": "b",
    "d": "h",
    "h": "d",
}


def reverse_complement(base_string):
    """
    Reverse complements the given base_string.
    :param base_string: String of bases to be reverse-complemented.
    :return: The reverse complement of the given base string.
    """

    return "".join(map(lambda b: RC_BASE_MAP[b], base_string[::-1]))


def get_confidence_factor(qual_string: str, scale_factor: float = CONF_FACTOR_SCALE) -> float:
    """Get the confidence factor for the given sequence to be tallied for use with STARCODE.
    quals are assumed to be phred scale quality scores in string format and will be converted to numerical values."""
    return scale_factor * reduce(
        operator.mul, map(lambda q: 1. - 10 ** (-(ord(q) - 33.) / 10), qual_string)
    )


def get_confidence_factor_raw_quals(quals: array.array, scale_factor: float = CONF_FACTOR_SCALE) -> float:
    """Get the confidence factor for the given sequence to be tallied for use with STARCODE.
    quals are assumed to be numerical already and will not be type converted."""
    return scale_factor * reduce(
        operator.mul, map(lambda q: 1. - 10 ** (-q/10), quals)
    )


def has_cbc(read):
    return read.has_tag(READ_RAW_BARCODE_TAG)


def has_umi(read):
    return read.has_tag(READ_RAW_UMI_TAG)


def get_model_name_from_bam_header(header):
    return get_model_name_from_bam_header(header)['name']


def get_model_from_bam_header(header):
    model_jsons = get_models_from_bam_header(header)

    if len(model_jsons) > 1:
        logger.warning(f"Loading model {model_jsons[0]['name']}, but more than one detected ({', '.join([m['name'] for m in model_jsons])}).")

    return model_jsons[0]


def get_models_from_bam_header(header):
    model_jsons = []
    for pg in header.as_dict()['PG']:
        try:
            if pg['PN'] == 'longbow' and (pg['ID'].startswith('longbow-annotate') or pg['ID'].startswith('longbow-pad')):
                desc, models_str = pg['DS'].split('MODEL(s): ')
                model_json = json.loads(models_str)
                model_jsons.append(model_json)
        except KeyError:
            continue

    return model_jsons


def bam_header_has_model(header):
    if 'PG' in header.as_dict():
        for pg in header.as_dict()['PG']:
            if pg['PN'] == 'longbow' and 'annotate' in pg['ID']:
                return True

    return False


def generate_read_name(movie_name, zmw, split_read_index):
    return f'{movie_name}/1{zmw:09d}{split_read_index:03d}/ccs'
