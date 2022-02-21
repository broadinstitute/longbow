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
def filtered_bam_file_from_pipeline(request):
    input_bam, model_name = request.param

    with tempfile.NamedTemporaryFile(delete=True) as annotate_bam, \
         tempfile.NamedTemporaryFile(delete=True) as filter_bam:

        runner = CliRunner()

        result_annotate = runner.invoke(longbow, ["annotate", "-m", model_name, "-f", "-o", annotate_bam.name, input_bam])
        assert result_annotate.exit_code == 0

        result_filter = runner.invoke(longbow, ["filter",   "-m", model_name, "-f", "-o", filter_bam.name,   annotate_bam.name])
        assert result_filter.exit_code == 0

        # Yield file here so that when we return, we get to clean up automatically
        yield filter_bam.name


def test_segment_from_file(tmpdir, filtered_bam_file_from_pipeline):
    actual_file = tmpdir.join(f"segment_actual_out.bam")
    args = ["segment", "-f", "-o", actual_file, filtered_bam_file_from_pipeline]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0


def test_segment_from_pipe(tmpdir, filtered_bam_file_from_pipeline):
    actual_file = tmpdir.join(f"filter_actual_out.pipe.bam")

    proc = subprocess.Popen(
        [ sys.executable, "-m", "longbow", "segment", "-f", "-o", actual_file ],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(filtered_bam_file_from_pipeline, proc)

    assert proc.returncode == 0
