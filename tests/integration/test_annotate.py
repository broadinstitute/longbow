import pytest

import subprocess
import sys
import os

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_reads_files_equal
from ..utils import cat_file_to_pipe

################################################################################

TOOL_NAME = "annotate"

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep
EXPECTED_DATA_FOLDER = TEST_DATA_FOLDER + TOOL_NAME + os.path.sep

################################################################################


@pytest.mark.parametrize("input_bam, expected_bam, model_name", [
    [TEST_DATA_FOLDER + "mas15_test_input.bam", EXPECTED_DATA_FOLDER + "mas15v2_expected.bam", "mas_15_sc_10x5p_single_none"],
    [TEST_DATA_FOLDER + "mas10_test_input.bam", EXPECTED_DATA_FOLDER + "mas10v2_expected.bam", "mas_10_sc_10x5p_single_none"],
])
def test_annotate(tmpdir, input_bam, expected_bam, model_name):

    actual_file = tmpdir.join(f"{TOOL_NAME}_actual_out.{model_name}.bam")
    args = ["annotate", "-t", 1, "-v", "INFO", "-m", model_name, input_bam, "-o", str(actual_file)]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    os.system(f"cp {str(actual_file)} .")

    assert result.exit_code == 0
    assert_reads_files_equal(actual_file, expected_bam, order_matters=True)


@pytest.mark.parametrize("input_bam, expected_bam, model_name", [
    [TEST_DATA_FOLDER + "mas15_test_input.bam", EXPECTED_DATA_FOLDER + "mas15v2_expected.bam", "mas_15_sc_10x5p_single_none"],
    [TEST_DATA_FOLDER + "mas10_test_input.bam", EXPECTED_DATA_FOLDER + "mas10v2_expected.bam", "mas_10_sc_10x5p_single_none"],
])
def test_annotate_from_pipe(tmpdir, input_bam, expected_bam, model_name):
    actual_bam = tmpdir.join(f"annotate_actual_out.{model_name}.pipe.bam")

    proc = subprocess.Popen(
        [sys.executable, "-m", "longbow", "annotate", "-m", model_name, "-f", "-o", actual_bam],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(input_bam, proc)

    assert proc.returncode == 0
    assert_reads_files_equal(actual_bam, expected_bam)
