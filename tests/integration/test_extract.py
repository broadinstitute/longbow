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
def segmented_bam_file_from_pipeline(request):
    input_bam, model_name = request.param

    with tempfile.NamedTemporaryFile(delete=True) as annotate_bam, \
         tempfile.NamedTemporaryFile(delete=True) as filter_bam, \
         tempfile.NamedTemporaryFile(delete=True) as segment_bam:

        runner = CliRunner()

        result_annotate = runner.invoke(longbow, ["annotate", "-t", 1, "-m", model_name, "-f", "-o", annotate_bam.name, str(input_bam)])
        assert result_annotate.exit_code == 0

        result_filter = runner.invoke(longbow, ["filter", "-m", model_name, "-f", "-o", filter_bam.name, annotate_bam.name])
        assert result_filter.exit_code == 0

        result_segment = runner.invoke(longbow, ["segment", "-t", 1, "-m", model_name, "-f", "-o", segment_bam.name,  filter_bam.name])
        assert result_segment.exit_code == 0

        # Yield file here so that when we return, we get to clean up automatically
        yield segment_bam.name


def test_extract_from_file(tmpdir, segmented_bam_file_from_pipeline):
    actual_file = tmpdir.join("extract_actual_out.bam")
    args = ["extract", "-f", "-o", actual_file, segmented_bam_file_from_pipeline]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0


def test_extract_from_pipe(tmpdir, segmented_bam_file_from_pipeline):
    actual_file = tmpdir.join("extract_actual_out.pipe.bam")

    args = ["extract", "-f", "-o", actual_file]

    runner = CliRunner()
    with open(segmented_bam_file_from_pipeline, "rb") as fh:
        result = runner.invoke(longbow, args, input=fh)

    assert result.exit_code == 0
