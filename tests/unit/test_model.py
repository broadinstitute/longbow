import pytest

from longbow.utils import model


@pytest.mark.parametrize("model_name, expected_value", [
    ["mas_15+sc_10x5p", True],
    ["mas_15+sc_10x3p", True],
    ["mas_15+bulk_10x5p", False],
    ["mas_10+sc_10x5p", True],
    ["mas_15+spatial_slideseq", False],
    ["mas_15+bulk_teloprimeV2", False],
    ["isoseq+sc_10x5p", True],
])
def test_has_cell_barcode_annotation(model_name, expected_value):
    lb_model = model.LibraryModel.build_pre_configured_model(model_name)
    assert lb_model.has_cell_barcode_annotation == expected_value


@pytest.mark.parametrize("model_name, expected_value", [
    ["mas_15+sc_10x5p", True],
    ["mas_15+sc_10x3p", True],
    ["mas_15+bulk_10x5p", True],
    ["mas_10+sc_10x5p", True],
    ["mas_15+spatial_slideseq", True],
    ["mas_15+bulk_teloprimeV2", False],
    ["isoseq+sc_10x5p", True],
])
def test_has_umi_annotation(model_name, expected_value):
    lb_model = model.LibraryModel.build_pre_configured_model(model_name)
    assert lb_model.has_umi_annotation == expected_value
