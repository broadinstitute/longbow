import pathlib

import pytest
from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data"


@pytest.mark.parametrize(
    "input_bam, model_exp",
    [
        [TEST_DATA_FOLDER / "mas15_test_input.bam", "mas_15+sc_10x5p"],
    ],
)
@pytest.mark.slow
def test_peek_from_file(tmpdir, input_bam, model_exp):
    actual_file = tmpdir.join(f"peek_actual_out.{model_exp}.txt")
    args = ["-t", 4, "peek", str(input_bam), "-f", "-o", actual_file]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0
    with open(actual_file, "r") as f:
        lines = f.readlines()
        assert lines[0].rstrip() == model_exp


@pytest.mark.parametrize(
    "input_bam, model_exp",
    [
        [TEST_DATA_FOLDER / "mas15_test_input.bam", "mas_15+sc_10x5p"],
    ],
)
@pytest.mark.slow
def test_peek_from_pipe(tmpdir, input_bam, model_exp):
    actual_file = tmpdir.join(f"peek_actual_out.{model_exp}.pipe.txt")
    args = ["-t", 4, "peek", "-f", "-o", str(actual_file)]

    runner = CliRunner()
    with open(input_bam, "rb") as fh:
        result = runner.invoke(longbow, args, input=fh)

    assert result.exit_code == 0
    with open(actual_file, "r") as f:
        lines = f.readlines()
        assert lines[0].rstrip() == model_exp
