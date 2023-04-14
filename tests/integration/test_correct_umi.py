import pathlib

import pytest
from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_reads_files_equal, convert_sam_to_bam

TOOL_NAME = "correct_umi"
TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data" / TOOL_NAME


@pytest.fixture(
    scope="function",
    params=[
        (
            TEST_DATA_FOLDER / "correct_umi_test_data.sam",
            TEST_DATA_FOLDER / "correct_umi_expected_corrected_data.sam",
            TEST_DATA_FOLDER / "correct_umi_expected_uncorrected_data.sam",
        )
    ],
)
def input_data_files(tmpdir, request):
    input_sam, expected_umi_corrected_sam, expected_umi_uncorrected_sam = request.param

    # Convert test files to bam:
    input_bam = tmpdir.join("input.bam")
    convert_sam_to_bam(input_sam, input_bam)

    expected_umi_corrected_bam = tmpdir.join("expected.bam")
    convert_sam_to_bam(expected_umi_corrected_sam, expected_umi_corrected_bam)

    return input_bam, expected_umi_corrected_bam, expected_umi_uncorrected_sam


def test_correct_umi(tmpdir, input_data_files):
    input_bam, expected_corrected_bam, expected_uncorrected_sam = input_data_files

    actual_umi_corrected_file = tmpdir.join(
        f"{TOOL_NAME}_actual_umi_corrected_out.mas15.bam"
    )
    actual_umi_uncorrected_file = tmpdir.join(
        f"{TOOL_NAME}_actual_umi_uncorrected_out.mas15.bam"
    )
    args = [
        "correct_umi",
        "-f",
        "--pre-extracted",
        "-l",
        10,
        "-o",
        str(actual_umi_corrected_file),
        "-x",
        str(actual_umi_uncorrected_file),
        str(input_bam),
    ]

    runner = CliRunner(mix_stderr=False)
    with runner.isolated_filesystem():
        result = runner.invoke(longbow, args)

    assert result.exit_code == 0

    # Equal files result as True:
    assert_reads_files_equal(
        actual_umi_corrected_file, expected_corrected_bam, order_matters=True
    )
    assert_reads_files_equal(
        actual_umi_uncorrected_file, expected_uncorrected_sam, order_matters=True
    )


def test_correct_umi_from_pipe(tmpdir, input_data_files):
    (
        input_bam,
        expected_umi_corrected_bam,
        expected_umi_uncorrected_sam,
    ) = input_data_files

    actual_umi_corrected_file = tmpdir.join(
        f"{TOOL_NAME}_actual_umi_corrected_out.mas15.pipe.bam"
    )
    actual_umi_uncorrected_file = tmpdir.join(
        f"{TOOL_NAME}_actual_umi_uncorrected_out.mas15.pipe.bam"
    )

    args = [
        "correct_umi",
        "-f",
        "--pre-extracted",
        "-l",
        10,
        "-o",
        str(actual_umi_corrected_file),
        "-x",
        str(actual_umi_uncorrected_file),
    ]

    runner = CliRunner()
    with runner.isolated_filesystem(), open(input_bam, "rb") as fh:
        result = runner.invoke(longbow, args, input=fh)

    assert result.exit_code == 0

    # Equal files result as True:
    assert_reads_files_equal(
        actual_umi_corrected_file, expected_umi_corrected_bam, order_matters=True
    )
    assert_reads_files_equal(
        actual_umi_uncorrected_file, expected_umi_uncorrected_sam, order_matters=True
    )
