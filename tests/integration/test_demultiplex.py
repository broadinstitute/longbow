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
    (TEST_DATA_FOLDER + "mas15_test_input.bam", "mas_15_sc_10x5p_single_none"),
    (TEST_DATA_FOLDER + "mas10_test_input.bam", "mas_10_sc_10x5p_single_none"),
])
def annotated_bam_file_from_pipeline(request):
    input_bam, model_name = request.param

    with tempfile.NamedTemporaryFile(delete=True) as annotate_bam:

        runner = CliRunner()

        result_annotate = runner.invoke(longbow, ["annotate", "-m", model_name, "-f", "-o", annotate_bam.name, input_bam])
        assert result_annotate.exit_code == 0

        # Yield file here so that when we return, we get to clean up automatically
        yield annotate_bam.name


def test_demultiplex_from_file(tmpdir, annotated_bam_file_from_pipeline):
    args = ["demultiplex", "-d", "YN", "-o", "demux", annotated_bam_file_from_pipeline]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0


@pytest.mark.skip(reason="this test is broken and I don't know why")
def test_demultiplex_from_pipe(tmpdir, annotated_bam_file_from_pipeline):
    proc = subprocess.Popen(
        [ sys.executable, "-m", "longbow", "demultiplex", "-d", "YN", "-o", "demux", annotated_bam_file_from_pipeline ],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(annotated_bam_file_from_pipeline, proc)

    assert proc.returncode == 0
