import pytest
import os
import sys
import subprocess

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import cat_file_to_pipe

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep


@pytest.mark.parametrize("input_bam, model_exp", [
    [TEST_DATA_FOLDER + "mas15_test_input.bam", "mas_15_sc_10x5p_single_none"],
])
def test_peek_from_file(tmpdir, input_bam, model_exp):

    actual_file = tmpdir.join(f"peek_actual_out.{model_exp}.txt")
    args = ["peek", "-t", 4, input_bam, "-f", "-o", actual_file]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0
    with open(actual_file, 'r') as f:
        l = f.readlines()
        assert l[0].rstrip() == model_exp


@pytest.mark.parametrize("input_bam, model_exp", [
    [TEST_DATA_FOLDER + "mas15_test_input.bam", "mas_15_sc_10x5p_single_none"],
])
def test_peek_from_pipe(tmpdir, input_bam, model_exp):
    actual_file = tmpdir.join(f"peek_actual_out.{model_exp}.pipe.txt")

    proc = subprocess.Popen(
        [ sys.executable, "-m", "longbow", "peek", "-f", "-o", actual_file ],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(input_bam, proc)

    assert proc.returncode == 0
    with open(actual_file, 'r') as f:
        l = f.readlines()
        assert l[0].rstrip() == model_exp
