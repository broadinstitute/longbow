import gzip
from construct import *


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