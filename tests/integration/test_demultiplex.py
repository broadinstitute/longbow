import pathlib

import pytest
from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data"
TEST_PARAMS = [
    TEST_DATA_FOLDER / "annotate" / "mas_15+sc_10x5p.expected.bam",
    TEST_DATA_FOLDER / "annotate" / "mas_10+sc_10x5p.expected.bam",
]


@pytest.mark.parametrize("input_bam", TEST_PARAMS)
def test_demultiplex_from_file(tmpdir, input_bam):
    args = ["demultiplex", "-d", "YN", "-o", "demux", str(input_bam)]

    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(longbow, args)

    assert result.exit_code == 0


@pytest.mark.parametrize("input_bam", TEST_PARAMS)
def test_demultiplex_from_pipe(tmpdir, input_bam):
    args = ["demultiplex", "-d", "YN", "-o" "demux"]

    runner = CliRunner()
    with runner.isolated_filesystem(), open(input_bam, "rb") as fh:
        result = runner.invoke(longbow, args, input=fh)

    assert result.exit_code == 0
