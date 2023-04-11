import pytest
import pathlib

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_reads_files_equal
from ..utils import convert_sam_to_bam


TOOL_NAME = "tagfix"
TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data" / TOOL_NAME
TEST_PARAMS = [
    [
        TEST_DATA_FOLDER / "tagfix_test_data.sam",
        TEST_DATA_FOLDER / "tagfix_test_data.expected.sam",
    ],
]


@pytest.mark.parametrize("input_sam, expected_sam", TEST_PARAMS)
def test_tagfix(tmpdir, input_sam, expected_sam):
    # Convert test files to bam:
    input_bam = tmpdir.join("input.bam")
    convert_sam_to_bam(input_sam, input_bam)

    expected_bam = tmpdir.join("expected.bam")
    convert_sam_to_bam(expected_sam, expected_bam)

    actual_bam_out = tmpdir.join(f"{TOOL_NAME}_actual_out.bam")
    args = [
        "-t",
        "1",
        "tagfix",
        "-v",
        "INFO",
        "-o",
        str(actual_bam_out),
        str(input_bam),
    ]

    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(longbow, args)
    assert result.exit_code == 0

    # Equal files result as True:
    assert_reads_files_equal(actual_bam_out, expected_bam, order_matters=True)


@pytest.mark.parametrize("input_sam, expected_sam", TEST_PARAMS)
def test_tagfix_from_pipe(tmpdir, input_sam, expected_sam):
    # Convert test files to bam:
    input_bam = tmpdir.join("input.bam")
    convert_sam_to_bam(input_sam, input_bam)

    expected_bam = tmpdir.join("expected.bam")
    convert_sam_to_bam(expected_sam, expected_bam)

    actual_bam_out = tmpdir.join(f"{TOOL_NAME}_actual_out.bam")

    args = ["-t", "1", "tagfix", "-v", "INFO", "-o", str(actual_bam_out)]

    runner = CliRunner()
    with open(input_bam, "rb") as fh:
        result = runner.invoke(longbow, args, input=fh)

    assert result.exit_code == 0

    # Equal files result as True:
    assert_reads_files_equal(actual_bam_out, expected_bam, order_matters=True)
