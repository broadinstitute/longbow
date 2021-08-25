import logging
import sys

import time

import click
import click_log
import tqdm

import pysam
import multiprocessing as mp

import numpy as np

from ..utils import bam_utils
from ..utils.model import LibraryModel

from ..annotate.command import SegmentInfo
from ..annotate.command import get_segments


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("stats")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-o",
    "--output-prefix",
    default="longbow_stats",
    type=str,
    help="prefix to give to output files",
)
@click.option(
    "-m",
    "--model",
    default="mas15",
    show_default=True,
    help="The model to use for annotation.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(threads, output_prefix, do_simple_splitting, model, input_bam):
    """Calculate and produce stats on the given input bam file."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        colour="green",
        file=sys.stderr,
        leave=False,
        disable=not sys.stdin.isatty(),
    ) as pbar:
        # Get our model:
        if LibraryModel.has_prebuilt_model(model):
            logger.info(f"Using %s", LibraryModel.pre_configured_models[model]["description"])
            m = LibraryModel.build_pre_configured_model(model)
        else:
            logger.info(f"Loading model from json file: %s", model)
            m = LibraryModel.from_json_file(model)

        # Create storage point for heatmap data:
        adapter_names = [array_element_adapters[0] for array_element_adapters in m.array_element_structure]
        ligation_heat_matrix = np.zeros((len(adapter_names) * 2, len(adapter_names) * 2), dtype=int)
        index_map = dict()
        for i, name in enumerate(adapter_names):
            index_map[name] = i
        for i, name in enumerate(adapter_names):
            index_map[name + "'"] = i + len(adapter_names)

        # Keep track of how many arrays start and end with each element:
        array_start_adapter_counts = {a: 0 for a in adapter_names}
        array_end_adapter_counts = {a: 0 for a in adapter_names}
        array_lengths = []

        for read in bam_file:
            # Get our read segments:
            try:
                _, segments = get_segments(read)
            except KeyError:
                logger.error(f"Input bam file does not contain longbow segmented reads!  "
                             f"No {bam_utils.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                sys.exit(1)

            print(segments)
            break

            # Update progress:
            pbar.update(1)

    logger.info("Processing statistics...")
    array_lengths = np.array(array_lengths)

    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)
