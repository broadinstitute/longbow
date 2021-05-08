import pytest
import os

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_bam_files_equal

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep
EXPECTED_DATA_FOLDER = TEST_DATA_FOLDER + "annotate" + os.path.sep


@pytest.mark.parametrize("input_bam, expected_bam, use_mas10", [
    [TEST_DATA_FOLDER + "mas15_test_input.bam", EXPECTED_DATA_FOLDER + "mas15_expected.bam", False],
    [TEST_DATA_FOLDER + "mas10_test_input.bam", EXPECTED_DATA_FOLDER + "mas10_expected.bam", True],
])
def test_annotate(tmpdir, input_bam, expected_bam, use_mas10):

    if use_mas10:
        actual_file = tmpdir.join("annotate_actual_out.mas10.bam")
        args = ["annotate", "-t", 4, "-v", "INFO", "--m10", input_bam, "-o", actual_file]
    else:
        actual_file = tmpdir.join("annotate_actual_out.mas15.bam")
        args = ["annotate", "-t", 4, "-v", "INFO", input_bam, "-o", actual_file]

    runner = CliRunner()
    runner.invoke(longbow, args)

    # Equal files result as True:
    assert_bam_files_equal(actual_file, expected_bam)


