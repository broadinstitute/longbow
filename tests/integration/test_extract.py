import pytest
import pathlib

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow


TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data"
TEST_PARAMS = [
    TEST_DATA_FOLDER / "segment" / "mas_15+sc_10x5p.expected.bam",
    TEST_DATA_FOLDER / "segment" / "mas_10+sc_10x5p.expected.bam",
]


@pytest.mark.parametrize("input_bam", TEST_PARAMS)
def test_extract_from_file(tmpdir, input_bam):
    actual_file = tmpdir.join("extract_actual_out.bam")
    args = ["extract", "-f", "-o", actual_file, str(input_bam)]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0


@pytest.mark.parametrize("input_bam", TEST_PARAMS)
def test_extract_from_pipe(tmpdir, input_bam):
    actual_file = tmpdir.join("extract_actual_out.pipe.bam")

    args = ["extract", "-f", "-o", actual_file]

    runner = CliRunner()
    with open(input_bam, "rb") as fh:
        result = runner.invoke(longbow, args, input=fh)

    assert result.exit_code == 0
