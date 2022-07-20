import pytest

import pysam
import tempfile
import os
import json

from longbow.utils import bam_utils
from longbow.utils import model


TEST_DATA_FOLDER = path = os.path.abspath(
    __file__ + os.path.sep + "../../" + os.path.sep + "test_data"
)


@pytest.mark.parametrize("model_name, expected_value", [
    ["mas_15_sc_10x5p_single_none", True],
    ["mas_15_sc_10x3p_single_none", True],
    ["mas_15_bulk_10x5p_single_internal", False],
    ["mas_10_sc_10x5p_single_none", True],
    ["mas_15_spatial_slide-seq_single_none", False],
    ["mas_15_bulk_teloprimeV2_single_none", False],
    ["isoseq_1_sc_10x5p_single_none", True],
])
def test_has_cell_barcode_annotation(model_name, expected_value):
    lb_model = model.LibraryModel.build_pre_configured_model(model_name)
    assert lb_model.has_cell_barcode_annotation == expected_value


@pytest.mark.parametrize("model_name, expected_value", [
    ["mas_15_sc_10x5p_single_none", True],
    ["mas_15_sc_10x3p_single_none", True],
    ["mas_15_bulk_10x5p_single_internal", True],
    ["mas_10_sc_10x5p_single_none", True],
    ["mas_15_spatial_slide-seq_single_none", True],
    ["mas_15_bulk_teloprimeV2_single_none", False],
    ["isoseq_1_sc_10x5p_single_none", True],
])
def test_has_umi_annotation(model_name, expected_value):
    lb_model = model.LibraryModel.build_pre_configured_model(model_name)
    assert lb_model.has_umi_annotation == expected_value
