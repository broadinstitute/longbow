import sys
import gzip
import inspect
import re
import logging
import click_log

from construct import *
from inspect import getframeinfo, currentframe, getdoc

import pysam

from ..meta import VERSION

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("bam_utils")
click_log.basic_config(logger)

PB_READ_NAME_RE = re.compile("m[0-9]+e?_[0-9]{6}_[0-9]{6}/[0-9]+/.*")

# Constants for bam file reading / writing:
SEGMENTS_TAG = "SG"
SEGMENTS_RC_TAG = "RC"
SEGMENT_TAG_DELIMITER = ","
READ_MODEL_NAME_TAG = "YN"
READ_MODEL_SCORE_TAG = "YS"
READ_IS_VALID_FOR_MODEL_TAG = "YV"
READ_FIRST_KEY_SEG_TAG = "YK"
READ_NUM_KEY_SEGMENTS_TAG = "YG"

READ_UMI_TAG = ""


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


def create_bam_header_with_program_group(command_name, base_bam_header, description=None):
    """Create a pysam.AlignmentHeader object with program group (PG) information populated by the given arguments.

    This function is intended to be called from the 'main' function of a longbow subcommand because it uses reflection
    to pull in the first line of the docstring from the main function as the description (DS field)."""

    bam_header_dict = base_bam_header.to_dict()

    if not description:
        prev_frame = currentframe().f_back
        description = getdoc(prev_frame.f_globals['main']).split("\n")[0]

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
        if os.path.exists(f):
            if exist_ok:
                logger.warning(f"Output file exists: {f}.  Overwriting.")
            else:
                logger.error(f"Output file already exists: {f}!")
                do_files_exist = True
    if do_files_exist:
        sys.exit(1)
