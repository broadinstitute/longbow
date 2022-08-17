import json
import logging
import time
import os
import sys

import click
import click_log

import tqdm
import pysam
from construct import *

import longbow.utils.constants
from ..utils import bam_utils
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel
from ..annotate.command import get_segments

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("sift")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-p",
    "--pbi",
    required=False,
    type=click.Path(exists=True),
    help="BAM .pbi index file",
)
@click.option(
    "-o",
    "--output-bam",
    default="-",
    type=click.Path(exists=False),
    help="filtered bam output (passing reads only)  [default: stdout]",
)
@click.option(
    "-x",
    "--reject-bam",
    default="/dev/null",
    type=click.Path(exists=False),
    help="Filtered bam output (failing reads only)  [default: /dev/null]",
)
@click.option(
    "-m",
    "--model",
    help="The model to use for annotation.  If not specified, it will be autodetected from "
         "the BAM header.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name "
         "and Longbow will attempt to read in the file and create a LibraryModel from it.  "
         "Longbow will assume the contents are the configuration of a LibraryModel as per "
         "LibraryModel.to_json()."
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, output_bam, reject_bam, model, force, input_bam):
    """Filter segmented reads by conformation to expected cDNA design."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("About to Filter %d reads", read_count)

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files([output_bam, reject_bam], exist_ok=force)

    logger.info(f"Writing reads that conform to the model to: {output_bam}")
    logger.info(f"Writing reads that do not conform to the model to: {reject_bam}")

    # Open our input bam file:
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file, \
            tqdm.tqdm(
            desc="Progress",
            unit=" read",
            colour="green",
            file=sys.stderr,
            disable=not sys.stdin.isatty(),
            total=read_count
            ) as pbar:

        # Get our model:
        if model is None:
            lb_model = LibraryModel.from_json_obj(bam_utils.get_model_from_bam_header(bam_file.header))
        elif model is not None and LibraryModel.has_prebuilt_model(model):
            lb_model = LibraryModel.build_pre_configured_model(model)
        else:
            lb_model = LibraryModel.from_json_file(model)

        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        # logger.info(f"Filtering according to {lb_model.name} model ordered key adapters: "
        #             f"{', '.join(lb_model.key_adapters)}")

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[lb_model])

        # Setup output files:
        with pysam.AlignmentFile(output_bam, "wb", header=out_header) as passing_bam_file, \
                pysam.AlignmentFile(reject_bam, "wb", header=out_header) as failing_bam_file:

            num_passed = 0
            num_failed = 0

            cached_models = {}

            for read in bam_file:
                # Get our read segments:
                try:
                    _, segments = get_segments(read)
                except KeyError:
                    logger.error(f"Input bam file does not contain longbow segmented reads!  "
                                 f"No {longbow.utils.constants.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                    sys.exit(1)

                # Annotate the read with the model that was used in its validation:
                read.set_tag(longbow.utils.constants.READ_MODEL_NAME_TAG, lb_model.name)

                # Create (or retrieve from cache) the forward-connected model of just this section of the read
                segment_names = [s.name for s in segments]
                sub_model_name = "model_" + "_".join(segment_names)

                if sub_model_name not in cached_models:
                    direct_connections = {}
                    for i in range(len(segment_names)-1):
                        direct_connections[segment_names[i]] = [segment_names[i+1]]

                    sub_model = LibraryModel.from_json_obj({
                        "name": "model_" + "_".join(segment_names),
                        "description": "model_" + "_".join(segment_names),
                        "version": "0.0.1",
                        "array_element_structure": [segment_names],
                        "adapters": {key: lb_model.adapter_dict[key] for key in segment_names},
                        "direct_connections" : direct_connections,
                        "start_element_names": [s.name for s in segments if s.name not in lb_model.named_random_segments],
                        "end_element_names": [s.name for s in segments if s.name not in lb_model.named_random_segments],
                        "named_random_segments": lb_model.named_random_segments,
                        "coding_region": lb_model.coding_region,
                        "annotation_segments": lb_model.annotation_segments
                    })

                    sub_model.build_forward_connected()

                    cached_models[sub_model_name] = sub_model

                sub_model = cached_models[sub_model_name]

                logp, ppath = sub_model.annotate(read.query_sequence)
                qpath = bam_utils.collapse_annotations(ppath)

                is_valid = len(segment_names) == len(qpath)

                read.set_tag(longbow.utils.constants.SEGMENTS_TAG, longbow.utils.constants.SEGMENT_TAG_DELIMITER.join([q.to_tag() for q in qpath]))
                read.set_tag(longbow.utils.constants.READ_IS_VALID_FOR_MODEL_TAG, is_valid)

                if is_valid:
                    logger.debug("Read is %s valid: %s: adapter pattern: %s",
                                 lb_model.name,
                                 read.query_name,
                                 read.get_tag(longbow.utils.constants.SEGMENTS_TAG))

                    passing_bam_file.write(read)
                    num_passed += 1
                else:
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug("Read is not %s valid: %s: adapter pattern: %s",
                                     lb_model.name,
                                     read.query_name,
                                     read.get_tag(longbow.utils.constants.SEGMENTS_TAG))

                    failing_bam_file.write(read)
                    num_failed += 1

                pbar.update(1)

    # Calc some stats:
    pct_reads_passing = 100 * num_passed / (num_passed + num_failed) if (num_passed + num_failed) > 0 else 0
    pct_reads_failing = 100 * num_failed / (num_passed + num_failed) if (num_passed + num_failed) > 0 else 0

    # Yell at the user:
    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)
    logger.info(f"Total Reads Processed: %d", num_passed + num_failed)
    logger.info(f"# Reads Passing Model Filter: %d (%2.2f%%)", num_passed, pct_reads_passing)
    logger.info(f"# Reads Failing Model Filter: %d (%2.2f%%)", num_failed, pct_reads_failing)
