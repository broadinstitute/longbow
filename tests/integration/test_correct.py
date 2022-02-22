import pytest
import os
import sys
import subprocess
import tempfile

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import cat_file_to_pipe

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep


@pytest.fixture(scope="module", params=[
    (TEST_DATA_FOLDER + "mas15_test_input.bam", "mas15v2"),
    (TEST_DATA_FOLDER + "mas10_test_input.bam", "mas10v2"),
])
def extracted_bam_file_from_pipeline(request):
    input_bam, model_name = request.param

    with tempfile.NamedTemporaryFile(delete=True) as annotate_bam, \
         tempfile.NamedTemporaryFile(delete=True) as filter_bam, \
         tempfile.NamedTemporaryFile(delete=True) as segment_bam, \
         tempfile.NamedTemporaryFile(delete=True) as extract_bam:

        runner = CliRunner()

        result_annotate = runner.invoke(longbow, ["annotate", "-m", model_name, "-f", "-o", annotate_bam.name, input_bam])
        assert result_annotate.exit_code == 0

        result_filter = runner.invoke(longbow, ["filter",   "-m", model_name, "-f", "-o", filter_bam.name,   annotate_bam.name])
        assert result_filter.exit_code == 0

        result_segment = runner.invoke(longbow, ["segment",  "-m", model_name, "-f", "-o", segment_bam.name,  filter_bam.name])
        assert result_segment.exit_code == 0

        result_extract = runner.invoke(longbow, ["extract",  "-m", model_name, "-f", "-o", extract_bam.name,  segment_bam.name])
        assert result_extract.exit_code == 0

        # Yield file here so that when we return, we get to clean up automatically
        yield extract_bam.name


def test_correct_from_file(tmpdir, extracted_bam_file_from_pipeline):
    actual_file = tmpdir.join(f"correct_actual_out.bam")
    args = ["correct", "-f", "-o", actual_file, extracted_bam_file_from_pipeline]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0


@pytest.mark.skip(reason="`correct` command currently does not accept data from a pipe")
def test_correct_from_pipe(tmpdir, extracted_bam_file_from_pipeline):
    actual_file = tmpdir.join(f"correct_actual_out.pipe.bam")

    proc = subprocess.Popen(
        [ sys.executable, "-m", "longbow", "correct", "-f", "-o", actual_file ],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(extracted_bam_file_from_pipeline, proc)

    assert proc.returncode == 0
