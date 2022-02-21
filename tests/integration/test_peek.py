import pytest
import os

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import assert_bam_files_equal

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep


@pytest.mark.parametrize("input_bam, model_exp", [
    [TEST_DATA_FOLDER + "mas15_test_input.bam", "mas15v2"],
])
def test_peek(tmpdir, input_bam, model_exp):

    actual_file = tmpdir.join("peek_actual_out.txt")
    args = ["peek", "-t", 4, input_bam, "-f", "-o", actual_file]

    runner = CliRunner()
    runner.invoke(longbow, args)

    with open(actual_file, 'r') as f:
        l = f.readlines()
        assert l[0].rstrip() == model_exp


