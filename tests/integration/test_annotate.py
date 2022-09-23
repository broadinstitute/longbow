import pytest
import pathlib
import os

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_reads_files_equal

################################################################################

TOOL_NAME = "annotate"

TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data"
EXPECTED_DATA_FOLDER = TEST_DATA_FOLDER / TOOL_NAME

################################################################################


@pytest.mark.parametrize("input_bam, expected_bam, model_name", [
    [TEST_DATA_FOLDER / "mas15_test_input.bam", EXPECTED_DATA_FOLDER / "mas15v2_expected.bam", "mas_15_sc_10x5p_single_none"],
    [TEST_DATA_FOLDER / "mas10_test_input.bam", EXPECTED_DATA_FOLDER / "mas10v2_expected.bam", "mas_10_sc_10x5p_single_none"],
])
def test_annotate(tmpdir, input_bam, expected_bam, model_name):

    actual_file = tmpdir.join(f"{TOOL_NAME}_actual_out.{model_name}.bam")
    args = ["annotate", "-t", 1, "-v", "INFO", "-m", model_name, str(input_bam), "-o", str(actual_file)]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    os.system(f"cp {str(actual_file)} .")

    assert result.exit_code == 0
    assert_reads_files_equal(actual_file, expected_bam, order_matters=True)


@pytest.mark.parametrize("input_bam, expected_bam, model_name", [
    [TEST_DATA_FOLDER / "mas15_test_input.bam", EXPECTED_DATA_FOLDER / "mas15v2_expected.bam", "mas_15_sc_10x5p_single_none"],
    [TEST_DATA_FOLDER / "mas10_test_input.bam", EXPECTED_DATA_FOLDER / "mas10v2_expected.bam", "mas_10_sc_10x5p_single_none"],
])
def test_annotate_from_pipe(tmpdir, input_bam, expected_bam, model_name):
    actual_bam = tmpdir.join(f"annotate_actual_out.{model_name}.pipe.bam")

    args = ["annotate", "-t", 1, "-v", "INFO", "-m", model_name, "-f", "-o", str(actual_bam), "-"]

    runner = CliRunner()
    with open(input_bam, "rb") as fh:
        result = runner.invoke(longbow, args, input=fh)

    assert result.exit_code == 0
    assert_reads_files_equal(actual_bam, expected_bam)
