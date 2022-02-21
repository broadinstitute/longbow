import pytest

import subprocess
import sys
import os
import time

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_bam_files_equal
from ..utils import cat_file_to_pipe

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep
EXPECTED_DATA_FOLDER = TEST_DATA_FOLDER + "annotate" + os.path.sep


@pytest.mark.parametrize("input_bam, expected_bam, model_name", [
    [TEST_DATA_FOLDER + "mas15_test_input.bam", EXPECTED_DATA_FOLDER + "mas15_expected.bam", "mas15v2"],
    [TEST_DATA_FOLDER + "mas10_test_input.bam", EXPECTED_DATA_FOLDER + "mas10_expected.bam", "mas10v2"],
])
def test_annotate_from_file(tmpdir, input_bam, expected_bam, model_name):
    actual_bam = tmpdir.join("annotate_actual_out.mas10.bam")
    args = ["annotate", "-t", 4, "-m", model_name, "-f", "-o", actual_bam, input_bam]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0
    assert_bam_files_equal(actual_bam, expected_bam)


@pytest.mark.parametrize("input_bam, expected_bam, model_name", [
    [TEST_DATA_FOLDER + "mas15_test_input.bam", EXPECTED_DATA_FOLDER + "mas15_expected.bam", "mas15v2"],
    [TEST_DATA_FOLDER + "mas10_test_input.bam", EXPECTED_DATA_FOLDER + "mas10_expected.bam", "mas10v2"],
])
def test_annotate_from_pipe(tmpdir, input_bam, expected_bam, model_name):
    actual_bam = tmpdir.join(f"annotate_actual_out.{model_name}.pipe.bam")

    proc = subprocess.Popen(
        [ sys.executable, "-m", "longbow", "annotate", "-m", model_name, "-f", "-o", actual_bam ],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(input_bam, proc)

    assert proc.returncode == 0
    assert_bam_files_equal(actual_bam, expected_bam)
