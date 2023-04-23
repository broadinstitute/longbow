from longbow.utils import model


def test_has_cell_barcode_annotation(builtin_model):
    lb_model = model.LibraryModel.build_pre_configured_model(builtin_model.name)
    assert lb_model.has_cell_barcode_annotation == builtin_model.has_bc


def test_has_umi_annotation(builtin_model):
    lb_model = model.LibraryModel.build_pre_configured_model(builtin_model.name)
    assert lb_model.has_umi_annotation == builtin_model.has_umi


def test_has_has_coding_region(builtin_model):
    lb_model = model.LibraryModel.build_pre_configured_model(builtin_model.name)
    assert lb_model.has_coding_region == builtin_model.has_coding


def test_has_named_random_segments(builtin_model):
    lb_model = model.LibraryModel.build_pre_configured_model(builtin_model.name)
    assert lb_model.has_named_random_segments == builtin_model.has_random


def test_num_array_elements(builtin_model):
    lb_model = model.LibraryModel.build_pre_configured_model(builtin_model.name)
    assert lb_model.num_array_elements == builtin_model.num_elem
