import logging
import os
import sys
import time

import click
import click_log
import pysam
import tqdm

import longbow.utils.constants

from ..utils import bam_utils
from ..utils.bam_utils import SegmentInfo
from ..utils.cli_utils import get_field_count_and_percent_string, zero_safe_div
from ..utils.constants import FFORMAT

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
    "LibraryModel.to_json().",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist.",
)
@click.argument(
    "input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb")
)
def main(pbi, output_bam, reject_bam, model, force, input_bam):
    """Filter reads by conformation to expected segment order."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
    if not read_count:
        read_count = bam_utils.get_read_count_from_bam_index(input_bam)

    if read_count:
        logger.info("About to Filter %d reads", read_count)

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files([output_bam, reject_bam], exist_ok=force)

    # Get our model:
    lb_model = bam_utils.load_model(model, input_bam)
    logger.info(f"Using {lb_model.name}: {lb_model.description}")

    logger.info(f"Writing reads that conform to the model to: {output_bam}")
    logger.info(f"Writing reads that do not conform to the model to: {reject_bam}")

    # Open our input bam file:
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file:

        logger.info(
            f"Filtering according to {lb_model.name} model ordered key adapters: "
            f"{', '.join(lb_model.key_adapters)}"
        )

        # Get our header from the input bam file:
        out_header = pysam.AlignmentHeader.from_dict(
            bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header)
        )

        # Setup output files:
        with pysam.AlignmentFile(
            output_bam, "wb", header=out_header
        ) as passing_bam_file, pysam.AlignmentFile(
            reject_bam, "wb", header=out_header
        ) as failing_bam_file:

            num_passed = 0
            num_failed = 0

            tot_num_valid_adapters = 0
            tot_num_failed_adapters = 0

            for read in tqdm.tqdm(
                bam_file,
                desc="Progress",
                unit=" read",
                colour="green",
                file=sys.stderr,
                disable=not sys.stdin.isatty(),
                total=read_count,
            ):
                # Get our read segments:
                try:
                    # Create SegmentInfo objects so we can deal with them better:
                    segments = tuple(
                        [
                            SegmentInfo.from_tag(s)
                            for s in read.get_tag(
                                longbow.utils.constants.SEGMENTS_TAG
                            ).split(longbow.utils.constants.SEGMENT_TAG_DELIMITER)
                        ]
                    )
                except KeyError:
                    logger.error(
                        f"Input bam file does not contain longbow segmented reads!  "
                        f"No {longbow.utils.constants.SEGMENTS_TAG} tag detected on read {read.query_name} !"
                    )
                    sys.exit(1)

                # Annotate the read with the model that was used in its validation:
                read.set_tag(longbow.utils.constants.READ_MODEL_NAME_TAG, lb_model.name)

                # Check to see if the read is valid by this model and write it out:
                (
                    is_valid,
                    num_valid_adapters,
                    first_valid_adapter_index,
                ) = lb_model.validate_segment_order(segments)

                if is_valid:
                    logger.debug(
                        "Read is %s valid: %s: first key adapter: [%d, %s], # key adapters: %d",
                        lb_model.name,
                        read.query_name,
                        first_valid_adapter_index,
                        lb_model.key_adapters[first_valid_adapter_index],
                        num_valid_adapters,
                    )

                    read.set_tag(
                        longbow.utils.constants.READ_IS_VALID_FOR_MODEL_TAG, True
                    )
                    read.set_tag(
                        longbow.utils.constants.READ_NUM_KEY_SEGMENTS_TAG,
                        num_valid_adapters,
                    )
                    read.set_tag(
                        longbow.utils.constants.READ_FIRST_KEY_SEG_TAG,
                        lb_model.key_adapters[first_valid_adapter_index],
                    )
                    passing_bam_file.write(read)
                    tot_num_valid_adapters += num_valid_adapters
                    num_passed += 1
                else:
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug(
                            "Read is not %s valid: %s: first key adapter: [%d, %s], # key adapters: %d, "
                            "key adapters detected: %s",
                            lb_model.name,
                            read.query_name,
                            first_valid_adapter_index,
                            lb_model.key_adapters[first_valid_adapter_index],
                            num_valid_adapters,
                            ",".join([s.name for s in segments]),
                        )

                    read.set_tag(
                        longbow.utils.constants.READ_IS_VALID_FOR_MODEL_TAG, False
                    )
                    failing_bam_file.write(read)
                    tot_num_failed_adapters += num_valid_adapters
                    num_failed += 1

    # Calc some stats:
    total_reads = num_passed + num_failed

    # Yell at the user:
    logger.info(f"Done. Elapsed time: %{FFORMAT}s.", time.time() - t_start)
    logger.info(f"Total Reads Processed: {num_passed + num_failed}")

    count_str, pct_str = get_field_count_and_percent_string(
        num_passed, total_reads, FFORMAT
    )
    logger.info(f"# Reads Passing Model Filter: {count_str} {pct_str}")

    count_str, pct_str = get_field_count_and_percent_string(
        num_failed, total_reads, FFORMAT
    )
    logger.info(f"# Reads Failing Model Filter: {count_str} {pct_str}")

    logger.info(
        f"Total # correctly ordered key adapters in passing reads: {tot_num_valid_adapters}"
    )
    logger.info(
        f"Total # correctly ordered key adapters in failing reads: {tot_num_failed_adapters}"
    )
    logger.info(
        f"Avg # correctly ordered key adapters per passing read: %{FFORMAT} [%d]",
        zero_safe_div(tot_num_valid_adapters, num_passed),
        len(lb_model.key_adapters),
    )
    logger.info(
        f"Avg # correctly ordered key adapters per failing read: %{FFORMAT} [%d]",
        zero_safe_div(tot_num_failed_adapters, num_passed),
        len(lb_model.key_adapters),
    )
