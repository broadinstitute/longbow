import logging
import pickle

import click
import matplotlib.pyplot as plt
import networkx as nx

from ..utils import cli_utils
from ..utils.model import MODEL_NAME_REGEX, LibraryModel
from ..utils.model_utils import ModelBuilder

logger = logging.getLogger(__name__)


@click.command("models")
@click.option(
    "-l",
    "--list-models",
    cls=cli_utils.MutuallyExclusiveOption,
    mutually_exclusive=["dump"],
    required=False,
    is_flag=True,
    default=False,
    help="List the names of all models supported natively by this version of Longbow.",
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
    "including: json representation, dot file representations, transmission matrix, state emission json file.",
)
def main(list_models, dump):
    """Get information about built-in Longbow models."""

    if list_models or (not list_models and dump is None):
        print("Longbow includes the following models:")

        print("\nArray models\n============")
        print_models(ModelBuilder.pre_configured_models["array"])

        print("\ncDNA models\n===========")
        print_models(ModelBuilder.pre_configured_models["cdna"])

        print(
            "\nSpecify a fully combined model via '<array model>+<cDNA model>' syntax, e.g. 'mas_15+sc_10x5p'."
        )

    if dump is not None:
        model_name = dump
        (array_model_name, cdna_model_name) = MODEL_NAME_REGEX.split(model_name, 2)

        if array_model_name not in ModelBuilder.pre_configured_models["array"].keys():
            logger.error(
                f"Given model name '{array_model_name}' not in preconfigured array models "
                f"{list(ModelBuilder.pre_configured_models['array'].keys())}"
            )

        if cdna_model_name not in ModelBuilder.pre_configured_models["cdna"].keys():
            logger.error(
                f"Given model name '{cdna_model_name}' not in preconfigured cDNA models "
                f"{list(ModelBuilder.pre_configured_models['cdna'].keys())}"
            )

        logger.info(f"Generating model: {model_name}")
        lb = LibraryModel(
            array_model=ModelBuilder.pre_configured_models["array"][array_model_name],
            cdna_model=ModelBuilder.pre_configured_models["cdna"][cdna_model_name],
            model_name=model_name,
        )

        logger.info(f"Dumping {lb.name}: {lb.description}")

        model_dump_base_name = (
            f"longbow_model-{lb.name}-Av{lb.array_version}_Cv{lb.cdna_version}"
        )

        logger.info(f"Dumping dotfile: {model_dump_base_name + '.dot'}")
        lb.dump_as_dotfile(
            out_file_name=model_dump_base_name + ".dot", do_subgraphs=False
        )

        logger.info(f"Dumping simple dotfile: {model_dump_base_name + '.simple.dot'}")
        lb.dump_as_dotfile_simple(out_file_name=model_dump_base_name + ".simple.dot")

        logger.info(
            f"Dumping json model specification: {model_dump_base_name + '.spec.json'}"
        )
        lb.to_json(f"{model_dump_base_name}.spec.json")

        logger.info(
            f"Dumping dense transition matrix: {model_dump_base_name + '.dense_transition_matrix.pickle'}"
        )
        with open(f"{model_dump_base_name}.dense_transition_matrix.pickle", "wb") as f:
            pickle.dump(lb.hmm.dense_transition_matrix(), f)

        logger.info(
            f"Dumping emission distributions: {model_dump_base_name + '.emission_distributions.txt'}"
        )
        with open(f"{model_dump_base_name}.emission_distributions.txt", "w") as f:
            print(lb.hmm, file=f, flush=True)

        logger.info(f"Creating model graph from {len(lb.hmm.states)} states...")
        labelsdict = {}
        for s in lb.hmm.states:
            if s.name.endswith("-start") or s.name.endswith("-end"):
                labelsdict[s] = s.name
            else:
                labelsdict[s] = ""

        layout = nx.spring_layout(lb.hmm.graph)

        logger.info("Rendering graph ...")
        nx.draw(lb.hmm.graph, with_labels=True, labels=labelsdict, pos=layout)

        logger.info(f"Writing model graph now to {model_dump_base_name}.graph.png ...")
        plt.savefig(f"{model_dump_base_name}.graph.png", dpi=300, bbox_inches="tight")


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
    print(
        f"{model_name:{table_col_widths[0] + 4}s}{version:{table_col_widths[1] + 4}s}{description}"
    )
    for model_name, version, description in model_info_list:
        print(
            f"{model_name:{table_col_widths[0]+4}s}{version:{table_col_widths[1]+4}s}{description}"
        )
