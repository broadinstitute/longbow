import pytest
from click.testing import CliRunner

import pathlib
import tempfile
import gzip

import pysam

from longbow.__main__ import main_entry as longbow


TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data"


@pytest.mark.parametrize(
    "input_bam",
    [
        TEST_DATA_FOLDER / "mas15_test_input.bam",
        TEST_DATA_FOLDER / "mas10_test_input.bam",
    ],
)
def test_convert_from_file(tmpdir, input_bam):
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(
        input_bam, "rb", require_index=False, check_sq=False
    ) as input_file:
        with tempfile.NamedTemporaryFile(delete=True, suffix=".fq.gz") as tf:
            with gzip.open(tf.name, "wb") as output_fq:
                for bam_record in input_file:
                    output_fq.write(
                        f"@{bam_record.query_name}\n{bam_record.query_sequence}\n+\n{bam_record.qual}\n".encode(
                            "utf-8"
                        )
                    )

            actual_bam = tmpdir.join("convert_actual_out.bam")
            args = ["convert", "-f", "-o", actual_bam, str(input_bam)]

            runner = CliRunner()
            result = runner.invoke(longbow, args)

            assert result.exit_code == 0
