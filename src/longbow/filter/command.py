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
    '--m10',
    is_flag=True,
    default=False,
    show_default=True,
    help="Use the 10 array element MAS-seq model."
)
@click.option(
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, out_prefix, m10, force, input_bam):
    """Filter reads by whether they conform to expected segment order."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    if m10:
        logger.info("Using MAS-seq 10 array element annotation model.")
        lb_model = LibraryModel.build_and_return_mas_seq_10_model()
        model_name = "mas10"
    else:
        logger.info("Using MAS-seq default annotation model.")
        lb_model = LibraryModel.build_and_return_mas_seq_model()
        model_name = "mas15"

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("About to Filter %d reads", read_count)

    # Create some output file names:
    passing_out_name = f"{out_prefix}_longbow_filter_passed.bam"
    failing_out_name = f"{out_prefix}_longbow_filter_failed.bam"

    # Check to see if the output files exist:
    for f in [passing_out_name, failing_out_name]:
        if os.path.exists(f):
            if force:
                logger.warning(f"Output file exists: {f}.  Overwriting.")
            else:
                logger.error(f"Output file exists: {f}.")
                sys.exit(1)

    logger.info(f"Writing reads that conform to the model to: {passing_out_name}")
    logger.info(f"Writing reads that do not conform to the model to: {failing_out_name}")

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
        out_header = bam_utils.create_bam_header_with_program_group("filter", bam_file.header)

        # Setup output files:
        with pysam.AlignmentFile(passing_out_name, "wb", header=out_header) as passing_bam_file, \
                pysam.AlignmentFile(failing_out_name, "wb", header=out_header) as failing_bam_file:

            for read in bam_file:
                # Get our read segments:
                try:
                    segments = get_segments(read)
                except KeyError:
                    logger.error(f"Input bam file does not contain longbow segmented reads!  "
                                 f"No {bam_utils.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                    sys.exit(1)

                # Annotate the read with the model that was used in its validation:
                read.set_tag(bam_utils.READ_MODEL_NAME_TAG, model_name)

                # Check to see if the read is valid by this model and write it out:
                if lb_model.validate_segment_order([s.name for s in segments]):
                    read.set_tag(bam_utils.READ_IS_VALID_FOR_MODEL_TAG, True)
                    passing_bam_file.write(read)
                else:
                    read.set_tag(bam_utils.READ_IS_VALID_FOR_MODEL_TAG, False)
                    failing_bam_file.write(read)

                pbar.update(1)

    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)
