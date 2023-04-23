import importlib.resources
from dataclasses import dataclass

import pytest

from longbow.utils import model_utils


@dataclass
class TestModel:
    name: str
    num_elem: int
    has_bc: bool
    has_umi: bool
    has_coding: bool
    has_random: bool


# all combinations of the models included in src/longbow/models
BUILTIN_MODELS = [
    TestModel(f"{a_name}+{c_name}", num_elem, *c_prop)
    for a_name, num_elem in [
        ("isoseq", 2),
        ("mas_15", 16),
        ("mas_16", 17),
        ("mas_10", 11),
    ]
    for c_name, c_prop in [
        ("bulk_10x5p", (False, True, True, True)),
        ("spatial_slideseq", (False, True, True, True)),
        ("sc_10x5p", (True, True, True, True)),
        ("sc_10x3p", (True, True, True, True)),
        ("bulk_teloprimeV2", (False, False, True, True)),
    ]
]


@pytest.fixture(scope="package", params=BUILTIN_MODELS)
def builtin_model(request):
    return request.param


@pytest.fixture(scope="package", autouse=True)
def preload_models():
    with importlib.resources.files("longbow.models") as model_dir:
        model_utils.load_models(model_dir)
