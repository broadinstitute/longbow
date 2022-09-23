import pytest
import pathlib
import tempfile

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data"


@pytest.fixture(scope="module", params=[
    (TEST_DATA_FOLDER / "mas15_test_input.bam", "mas_15_sc_10x5p_single_none"),
    (TEST_DATA_FOLDER / "mas10_test_input.bam", "mas_10_sc_10x5p_single_none"),
])
def annotated_bam_file_from_pipeline(request):
    input_bam, model_name = request.param

    with tempfile.NamedTemporaryFile(delete=True) as annotate_bam:

        runner = CliRunner()

        result_annotate = runner.invoke(longbow, ["annotate", "-m", model_name, "-f", "-o", annotate_bam.name, str(input_bam)])
        assert result_annotate.exit_code == 0

        # Yield file here so that when we return, we get to clean up automatically
        yield annotate_bam.name


def test_demultiplex_from_file(tmpdir, annotated_bam_file_from_pipeline):
    args = ["demultiplex", "-d", "YN", "-o", "demux", annotated_bam_file_from_pipeline]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0


def test_demultiplex_from_pipe(tmpdir, annotated_bam_file_from_pipeline):
    args = ["demultiplex", "-d", "YN", "-o" "demux"]

    runner = CliRunner()
    with open(annotated_bam_file_from_pipeline, "rb") as fh:
        result = runner.invoke(longbow, args, input=fh)

    assert result.exit_code == 0
