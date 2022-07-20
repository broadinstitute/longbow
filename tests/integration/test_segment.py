import pytest
import os
import sys
import subprocess
import tempfile

from click.testing import CliRunner

from longbow.__main__ import main_entry as longbow

from ..utils import cat_file_to_pipe
from ..utils import convert_sam_to_bam
from ..utils import assert_reads_files_equal

################################################################################

TOOL_NAME = "segment"

TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
) + os.path.sep
EXPECTED_DATA_FOLDER = TEST_DATA_FOLDER + TOOL_NAME + os.path.sep

################################################################################


@pytest.fixture(scope="module", params=[
    (TEST_DATA_FOLDER + "mas15_test_input.bam",
     EXPECTED_DATA_FOLDER + "mas_15_sc_10x5p_single_none.expected.bam",
     "mas_15_sc_10x5p_single_none"),

    (TEST_DATA_FOLDER + "mas10_test_input.bam",
     EXPECTED_DATA_FOLDER + "mas_10_sc_10x5p_single_none.expected.bam",
     "mas_10_sc_10x5p_single_none"),

    (TEST_DATA_FOLDER + "mas_15_bulk_10x5p_single_internal.bam",
     EXPECTED_DATA_FOLDER + "mas_15_bulk_10x5p_single_internal.expected.bam",
     "mas_15_bulk_10x5p_single_internal"),
])
def filtered_bam_file_from_pipeline(request):
    input_file, expected_file, model_name = request.param

    if input_file.endswith(".bam"):
        input_bam = input_file
    else:
        # Convert test file to bam:
        with tempfile.NamedTemporaryFile(delete=True) as input_bam:
            convert_sam_to_bam(input_file, input_bam)

    with tempfile.NamedTemporaryFile(delete=True) as annotate_bam, \
         tempfile.NamedTemporaryFile(delete=True) as filter_bam:

        runner = CliRunner()

        result_annotate = runner.invoke(longbow, ["annotate", "-t", 1, "-m", model_name,
                                                  "-f", "-o", annotate_bam.name, input_bam])
        assert result_annotate.exit_code == 0

        result_filter = runner.invoke(longbow, ["filter", "-m", model_name,
                                                "-f", "-o", filter_bam.name, annotate_bam.name])
        assert result_filter.exit_code == 0

        # Yield file here so that when we return, we get to clean up automatically
        yield filter_bam.name, expected_file, model_name


def test_segment_from_file(tmpdir, filtered_bam_file_from_pipeline):
    input_bam_file, expected_bam, model_name = filtered_bam_file_from_pipeline

    actual_file = tmpdir.join(f"segment_actual_out.bam")
    args = ["segment", "-t", 1, "-f", "-o", actual_file, input_bam_file]

    runner = CliRunner()
    result = runner.invoke(longbow, args)

    assert result.exit_code == 0
    assert_reads_files_equal(actual_file, expected_bam, order_matters=True)


def test_segment_from_pipe(tmpdir, filtered_bam_file_from_pipeline):
    input_bam_file, expected_bam, model_name = filtered_bam_file_from_pipeline

    actual_file = tmpdir.join(f"filter_actual_out.pipe.bam")

    proc = subprocess.Popen(
        [sys.executable, "-m", "longbow", "segment", "-t", "1", "-f", "-o", actual_file],
        stdin=subprocess.PIPE
    )

    cat_file_to_pipe(input_bam_file, proc)

    assert proc.returncode == 0
    assert_reads_files_equal(actual_file, expected_bam, order_matters=True)
