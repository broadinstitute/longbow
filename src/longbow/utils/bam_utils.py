import sys
import gzip
from construct import *
from inspect import getframeinfo, currentframe, getdoc

import pysam

from ..meta import VERSION


# Constants for bam file reading / writing:
SEGMENTS_TAG = "SG"
SEGMENTS_RC_TAG = "RC"
SEGMENT_TAG_DELIMITER = ","
READ_MODEL_NAME_TAG = "YN"
READ_MODEL_SCORE_TAG = "YS"
READ_IS_VALID_FOR_MODEL_TAG = "YV"


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

    This function is intended to be called from the main function of a longbow subcommand because it uses reflection to
    pull in the first line of the docstring from the main method as the description (DS field)."""

    bam_header_dict = base_bam_header.to_dict()

    # Add our program group to it:
    pg_dict = {
        "ID": f"longbow-{command_name}-{VERSION}",
        "PN": "longbow",
        "VN": f"{VERSION}",
        # Use reflection to get the first line of the doc string the caller - the main function for our header:
        "DS": description if description else getdoc(globals()[getframeinfo(currentframe()).function]).split("\n")[0],
        "CL": " ".join(sys.argv),
    }
    if "PG" in bam_header_dict:
        bam_header_dict["PG"].append(pg_dict)
    else:
        bam_header_dict["PG"] = [pg_dict]
    out_header = pysam.AlignmentHeader.from_dict(bam_header_dict)

    return out_header
