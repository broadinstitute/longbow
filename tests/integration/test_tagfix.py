import pytest
import os
import sys
import subprocess

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_reads_files_equal
from ..utils import cat_file_to_pipe
from ..utils import convert_sam_to_bam

################################################################################

TOOL_NAME = "tagfix"

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep + TOOL_NAME + os.path.sep


################################################################################


@pytest.mark.parametrize("input_sam, expected_sam", [
    [
        TEST_DATA_FOLDER + "tagfix_test_data.sam",
        TEST_DATA_FOLDER + "tagfix_test_data.expected.sam"
    ],
])
def test_tagfix(tmpdir, input_sam, expected_sam):

    # Convert test files to bam:
    input_bam = tmpdir.join("input.bam")
    convert_sam_to_bam(input_sam, input_bam)

    expected_bam = tmpdir.join("expected.bam")
    convert_sam_to_bam(expected_sam, expected_bam)

    actual_bam_out = tmpdir.join(f"{TOOL_NAME}_actual_out.bam")
    args = ["tagfix", "-t", 1, "-v", "INFO", "-o", str(actual_bam_out), str(input_bam)]

    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(longbow, args)
    assert result.exit_code == 0

    # Equal files result as True:
    assert_reads_files_equal(actual_bam_out, expected_bam, order_matters=True)


@pytest.mark.skip(reason="tagfix pipe input test is not working properly and needs to be debugged.")
@pytest.mark.parametrize("input_sam, expected_sam", [
    [
        TEST_DATA_FOLDER + "tagfix_test_data.sam",
        TEST_DATA_FOLDER + "tagfix_test_data.expected.sam"
    ],
])
def test_tagfix_from_pipe(tmpdir, input_sam, expected_sam):

    # Convert test files to bam:
    input_bam = tmpdir.join("input.bam")
    convert_sam_to_bam(input_sam, input_bam)

    expected_bam = tmpdir.join("expected.bam")
    convert_sam_to_bam(expected_sam, expected_bam)

    actual_bam_out = tmpdir.join(f"{TOOL_NAME}_actual_out.bam")

    proc = subprocess.Popen(
        [sys.executable, "-m", "longbow", "tagfix", "-t", 1, "-v", "INFO", "-o", str(actual_bam_out)],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(str(input_bam), proc)

    assert proc.returncode == 0

    # Equal files result as True:
    assert_reads_files_equal(actual_bam_out, expected_bam, order_matters=True)
