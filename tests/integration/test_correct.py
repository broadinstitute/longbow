import pytest
import os
import sys
import subprocess

import pysam

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_bam_files_equal
from ..utils import cat_file_to_pipe

################################################################################

TOOL_NAME = "correct"

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep + TOOL_NAME + os.path.sep

################################################################################


def convert_sam_to_bam(sam_path, out_bam_path):
    with pysam.AlignmentFile(sam_path, "r", check_sq=False, require_index=False) as input_file:
        with pysam.AlignmentFile(out_bam_path, "wb", header=input_file.header) as out_bam_file:
            for r in input_file:
                out_bam_file.write(r)

################################################################################


@pytest.mark.parametrize("input_sam, expected_bc_corrected_sam, expected_bc_uncorrected_sam", [
    [TEST_DATA_FOLDER + "correct_test_data.sam", TEST_DATA_FOLDER + "correct_expected_corrected_data.sam",
     TEST_DATA_FOLDER + "correct_expected_uncorrected_data.sam"],
])
def test_correct(tmpdir, input_sam, expected_bc_corrected_sam, expected_bc_uncorrected_sam):

    # Convert test files to bam:
    input_bam = tmpdir.join("input.bam")
    convert_sam_to_bam(input_sam, input_bam)

    expected_bc_corrected_bam = tmpdir.join("expected.bam")
    convert_sam_to_bam(expected_bc_corrected_sam, expected_bc_corrected_bam)

    actual_bc_corrected_file = tmpdir.join(f"{TOOL_NAME}_actual_out.mas15.bam")
    actual_bc_uncorrected_file = tmpdir.join(f"{TOOL_NAME}_actual_bc_uncorrected_out.mas15.bam")
    args = ["correct", "-t", 1, "-m", "mas15v2", "-v", "INFO",
            "-a", f"{TEST_DATA_FOLDER}barcode_allow_list.txt", str(input_bam), "-o", str(actual_bc_corrected_file),
            "--barcode-uncorrectable-bam", str(actual_bc_uncorrected_file)]

    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0

    os.system(f"cp -v {str(actual_bc_corrected_file)} .")

    # Equal files result as True:
    assert_bam_files_equal(actual_bc_corrected_file, expected_bc_corrected_bam, order_matters=True)
    assert_bam_files_equal(actual_bc_uncorrected_file, expected_bc_uncorrected_sam, order_matters=True)


@pytest.mark.skip(reason="`correct` command currently does not accept data from a pipe")
def test_correct_from_pipe(tmpdir, extracted_bam_file_from_pipeline):
    actual_file = tmpdir.join(f"correct_actual_out.pipe.bam")

    proc = subprocess.Popen(
        [sys.executable, "-m", "longbow", "correct", "-f", "-o", actual_file],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(extracted_bam_file_from_pipeline, proc)

    assert proc.returncode == 0
