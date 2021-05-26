import logging
import time
import os
import sys

import click
import click_log

import tqdm
import pysam
from construct import *

from ..utils import bam_utils
from ..utils.model import LibraryModel
from ..annotate.command import get_segments

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("filter")
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
    "--out-prefix",
    default=".",
    required=True,
    type=str,
    help="Output file prefix",
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
@click.option(
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, out_prefix, model, force, input_bam):
    """Filter reads by whether they conform to expected segment order."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Get our model:
    if LibraryModel.has_prebuilt_model(model):
        logger.info(f"Using %s", LibraryModel.pre_configured_models[model]["description"])
        lb_model = LibraryModel.build_pre_configured_model(model)
    else:
        logger.info(f"Loading model from json file: %s", model)
        lb_model = LibraryModel.from_json_file(model)

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("About to Filter %d reads", read_count)

    # Create some output file names:
    passing_out_name = f"{out_prefix}_longbow_filter_passed.bam"
    failing_out_name = f"{out_prefix}_longbow_filter_failed.bam"

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files([passing_out_name, failing_out_name], exist_ok=force)

    logger.info(f"Writing reads that conform to the model to: {passing_out_name}")
    logger.info(f"Writing reads that do not conform to the model to: {failing_out_name}")

    logger.info(f"Filtering according to {lb_model.name} model ordered key adapters: "
                f"{', '.join(lb_model.key_adapters)}")

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

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[lb_model])

        # Setup output files:
        with pysam.AlignmentFile(passing_out_name, "wb", header=out_header) as passing_bam_file, \
                pysam.AlignmentFile(failing_out_name, "wb", header=out_header) as failing_bam_file:

            num_passed = 0
            num_failed = 0

            tot_num_valid_adapters = 0
            tot_num_failed_adapters = 0

            for read in bam_file:
                # Get our read segments:
                try:
                    _, segments = get_segments(read)
                except KeyError:
                    logger.error(f"Input bam file does not contain longbow segmented reads!  "
                                 f"No {bam_utils.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                    sys.exit(1)

                # Annotate the read with the model that was used in its validation:
                read.set_tag(bam_utils.READ_MODEL_NAME_TAG, lb_model.name)

                # Check to see if the read is valid by this model and write it out:
                segment_names = [s.name for s in segments]
                is_valid, num_valid_adapters, first_valid_adapter_index = \
                    lb_model.validate_segment_order(segment_names)

                if is_valid:
                    logger.debug("Read is %s valid: %s: first key adapter: [%d, %s], # key adapters: %d",
                                 lb_model.name,
                                 read.query_name,
                                 first_valid_adapter_index,
                                 lb_model.key_adapters[first_valid_adapter_index],
                                 num_valid_adapters)

                    read.set_tag(bam_utils.READ_IS_VALID_FOR_MODEL_TAG, True)
                    read.set_tag(bam_utils.READ_NUM_KEY_SEGMENTS_TAG, num_valid_adapters)
                    read.set_tag(bam_utils.READ_FIRST_KEY_SEG_TAG, lb_model.key_adapters[first_valid_adapter_index])
                    passing_bam_file.write(read)
                    tot_num_valid_adapters += num_valid_adapters
                    num_passed += 1
                else:
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug("Read is not %s valid: %s: first key adapter: [%d, %s], # key adapters: %d, "
                                     "key adapters detected: %s",
                                     lb_model.name,
                                     read.query_name,
                                     first_valid_adapter_index,
                                     lb_model.key_adapters[first_valid_adapter_index],
                                     num_valid_adapters,
                                     lb_model.extract_key_segment_names(segment_names))

                    read.set_tag(bam_utils.READ_IS_VALID_FOR_MODEL_TAG, False)
                    failing_bam_file.write(read)
                    tot_num_failed_adapters += num_valid_adapters
                    num_failed += 1

                pbar.update(1)

    # Calc some stats:
    pct_reads_passing = 100 * num_passed / (num_passed + num_failed) if (num_passed + num_failed) > 0 else 0
    pct_reads_failing = 100 * num_failed / (num_passed + num_failed) if (num_passed + num_failed) > 0 else 0
    mean_adapters_per_passing_read = tot_num_valid_adapters/num_passed if num_passed > 0 else 0
    mean_adapters_per_failing_read = tot_num_failed_adapters/num_failed if num_failed > 0 else 0

    # Yell at the user:
    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)
    logger.info(f"Total Reads Processed: %d", num_passed + num_failed)
    logger.info(f"# Reads Passing Model Filter: %d (%2.2f%%)", num_passed, pct_reads_passing)
    logger.info(f"# Reads Failing Model Filter: %d (%2.2f%%)", num_failed, pct_reads_failing)
    logger.info(f"Total # correctly ordered key adapters in passing reads: %d", tot_num_valid_adapters)
    logger.info(f"Total # correctly ordered key adapters in failing reads: %d", tot_num_failed_adapters)
    logger.info(f"Avg # correctly ordered key adapters per passing read: %2.4f [%d]",
                mean_adapters_per_passing_read, len(lb_model.key_adapters))
    logger.info(f"Avg # correctly ordered key adapters per failing read: %2.4f [%d]",
                mean_adapters_per_failing_read, len(lb_model.key_adapters))
