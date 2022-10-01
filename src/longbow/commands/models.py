import logging
import time
import sys
import re

import pickle

import click
import click_log

import multiprocessing as mp

import networkx as nx
from networkx.drawing.nx_pydot import write_dot

import matplotlib.pyplot as plt

# from ..utils import model as LongbowModel
from ..utils.model import LibraryModel
from ..utils.model_utils import ModelBuilder

from ..utils import cli_utils

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("models")
click_log.basic_config(logger)


@click.command(name="models")
@click_log.simple_verbosity_option(logger)
@click.option(
    "-l",
    "--list-models",
    cls=cli_utils.MutuallyExclusiveOption,
    mutually_exclusive=["dump"],
    required=False,
    is_flag=True,
    default=False,
    help="List the names of all models supported natively by this version of Longbow."
)
@click.option(
    "-d",
    "--dump",
    cls=cli_utils.MutuallyExclusiveOption,
    mutually_exclusive=["list_models"],
    required=False,
    type=str,
    help="Dump the details of a given model.  "
         "This command creates a set of files corresponding to the given model "
         "including: json representation, dot file representations, transmission matrix, state emission json file."
)
def main(list_models, dump):
    """Get information about built-in Longbow models."""

    # logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    if list_models or (not list_models and dump is None):
        print("Longbow includes the following models:")

        print("\nArray models\n============")
        print_models(ModelBuilder.pre_configured_array_models)

        print("\ncDNA models\n===========")
        print_models(ModelBuilder.pre_configured_cdna_models)

        print("\nSpecify a fully combined model via '<array model>+<cDNA model>' syntax, e.g. 'mas_15+sc_10x5p'.")

    if dump is not None:
        model_name = dump
        (array_model_name, cdna_model_name) = re.split('\+', model_name, 2)

        if array_model_name not in ModelBuilder.pre_configured_array_models.keys():
            logger.error(f"Given model name '{array_model_name}' not in preconfigured array models {list(ModelBuilder.pre_configured_array_models.keys())}")

        if cdna_model_name not in ModelBuilder.pre_configured_cdna_models.keys():
            logger.error(f"Given model name '{cdna_model_name}' not in preconfigured cDNA models {list(ModelBuilder.pre_configured_cdna_models.keys())}")

        lb = LibraryModel(
            array_model=ModelBuilder.pre_configured_array_models[array_model_name], 
            cdna_model=ModelBuilder.pre_configured_cdna_models[cdna_model_name]
        )

        labelsdict = {}
        for s in lb.hmm.states:
            labelsdict[s] = s.name

        # layout = nx.spring_layout(lb.hmm.graph)
        # nx.draw(lb.hmm.graph, with_labels=True, labels=labelsdict, pos=layout)
        # plt.savefig('plotgraph.png', dpi=300, bbox_inches='tight')
        # plt.show()

        # Get out model:
        # if LibraryModel.has_prebuilt_model(model_name):
        #     m = LibraryModel.build_pre_configured_model(model_name)
        # else:
        #     logger.info(f"Loading model from json file: %s", model_name)
        #     m = LibraryModel.from_json_file(model_name)

        # m = LibraryModel.build_pre_configured_hierarchical_model(model_name)

        # logger.info(f"Dumping %s: %s", model_name, m.description)

        # model_dump_base_name = f"longbow_model_{m.name}.v{m.version}"

        # logger.info(f"Dumping dotfile: {model_dump_base_name + '.dot'}")
        # m.dump_as_dotfile(out_file_name=model_dump_base_name + ".dot", do_subgraphs=False)

        # logger.info(f"Dumping simple dotfile: {model_dump_base_name + '.simple.dot'}")
        # m.dump_as_dotfile_simple(out_file_name=model_dump_base_name + ".simple.dot")

        # logger.info(f"Dumping json model specification: {model_dump_base_name + '.spec.json'}")
        # m.to_json(f"{model_dump_base_name}.spec.json")

        # logger.info(f"Dumping dense transition matrix: {model_dump_base_name + '.dense_transition_matrix.pickle'}")
        # with open(f"{model_dump_base_name}.dense_transition_matrix.pickle", 'wb') as f:
        #     pickle.dump(m.hmm.dense_transition_matrix(), f)

        # logger.info(f"Dumping emission distributions: {model_dump_base_name + '.emission_distributions.txt'}")
        # with open(f"{model_dump_base_name}.emission_distributions.txt", 'w') as f:
        #     print(m.hmm, file=f, flush=True)


def print_models(models):
    model_info_list = []
    table_col_widths = [0, 0, 0]
    for k, v in models.items():
        model_info_list.append((k, v["version"], v["description"]))

        # Make sure we can print things nice and tidy:
        if len(k) > table_col_widths[0]:
            table_col_widths[0] = len(k)
        if len(v["version"]) > table_col_widths[1]:
            table_col_widths[1] = len(v["version"])

    model_name = "Name"
    version = "Version"
    description = "Description"
    print(f"{model_name:{table_col_widths[0] + 4}s}{version:{table_col_widths[1] + 4}s}{description}")
    for model_name, version, description in model_info_list:
        print(f"{model_name:{table_col_widths[0]+4}s}{version:{table_col_widths[1]+4}s}{description}")
