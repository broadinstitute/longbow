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
from collections import Counter

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
    "-l",
    "--validation-model",
    default="10x_sc_10x5p_single_none",
    help="The model to use for cDNA validation."
)
@click.option(
    '-f',
    '--force',
    is_flag=True,
    default=False,
    show_default=True,
    help="Force overwrite of the output files if they exist."
)
@click.option(
    "-s",
    "--stats",
    default="/dev/null",
    type=click.Path(exists=False),
    help="Table describing the ways in which the reads do not conform to expectation (failing reads only)"
         "[default: /dev/null]",
)
@click.option(
    "-u",
    "--summary-stats",
    default="/dev/null",
    type=click.Path(exists=False),
    help="Table containing summary statistics for the sifted reads that are output."
         "[default: /dev/null]",
)
@click.option(
    "-k",
    "--ignore-list",
    required=False,
    type=click.Path(exists=True),
    help="Txt file containing a list of read names to ignore.",
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, output_bam, reject_bam, model, validation_model, force, stats, summary_stats, ignore_list, input_bam):
    """Filter segmented reads by conformation to expected cDNA design."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("About to Filter %d reads", read_count)

    reads_to_ignore = set()
    if ignore_list and os.path.exists(ignore_list):
        logger.info(f"Ingesting read ignore list: {ignore_list}")
        with open(ignore_list, 'r') as f:
            for line in f:
                read_name = line.strip()
                if len(read_name) > 0:
                    reads_to_ignore.add(read_name)
        logger.info(f"Num reads to ignore: {len(reads_to_ignore)}")

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

        # TODO: Currently only one model works.  FIX THIS.
        # Make sure we've specified the model that works:
        if validation_model != "10x_sc_10x5p_single_none":
            logger.error(f"Given validation model is not yet implemented with `longbow sift`: {validation_model}")
            sys.exit(1)

        if lb_model.name != "mas_15_sc_10x5p_single_none":
            logger.error(f"Given model is not yet implemented with `longbow sift`: {lb_model.name}")
            sys.exit(1)

        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        # Get our header from the input bam file:
        out_header = bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header, models=[lb_model])

        # Setup output files:
        with pysam.AlignmentFile(output_bam, "wb", header=out_header) as passing_bam_file, \
                pysam.AlignmentFile(reject_bam, "wb", header=out_header) as failing_bam_file, \
                open(stats, 'w') as stats_file:

            num_passed = 0
            num_failed = 0
            num_ignored = 0

            # Create the forward-connected model of just this section of the read
            sub_model = LibraryModel.build_pre_configured_model(validation_model)
            sub_model.build_forward_connected()

            stats_file.write('\t'.join(['read_name', 'rq', '5p_Adapter', 'CBC', 'UMI', 'SLS', 'cDNA', 'Poly_A', '3p_Adapter', 'SG']) + '\n')

            for read in bam_file:

                if read.query_name in reads_to_ignore:
                    logger.debug(f"Ignoring read: {read.query_name}")
                    num_ignored += 1
                    continue

                # Get our read segments:
                try:
                    _, segments = get_segments(read)
                except KeyError:
                    logger.error(f"Input bam file does not contain longbow segmented reads!  "
                                 f"No {longbow.utils.constants.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                    sys.exit(1)

                # Annotate the read with the model that was used in its validation:
                logp, ppath = sub_model.annotate(read.query_sequence)
                qpath = bam_utils.collapse_annotations(ppath)

                counts = Counter([q.to_tag().split(":")[0] for q in qpath])

                is_valid = (counts['5p_Adapter'] == 1 and
                            counts['CBC'] == 1 and
                            counts['UMI'] == 1 and
                            counts['SLS'] == 1 and
                            counts['cDNA'] == 1 and
                            counts['Poly_A'] == 1 and
                            counts['3p_Adapter'] == 1)

                read.set_tag(longbow.utils.constants.READ_MODEL_NAME_TAG, lb_model.name)
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

                    stats_file.write('\t'.join([
                        read.query_name,
                        str(read.get_tag("rq")),
                        str(counts['5p_Adapter']),
                        str(counts['CBC']),
                        str(counts['UMI']),
                        str(counts['SLS']),
                        str(counts['cDNA']),
                        str(counts['Poly_A']),
                        str(counts['3p_Adapter']),
                        longbow.utils.constants.SEGMENT_TAG_DELIMITER.join([q.to_tag() for q in qpath])
                        ]) + '\n'
                    )

                pbar.update(1)

    # Calc some stats:
    tot_reads = num_passed + num_failed + num_ignored
    pct_reads_passing = 100 * num_passed / tot_reads if tot_reads > 0 else 0
    pct_reads_failing = 100 * num_failed / tot_reads if tot_reads > 0 else 0
    pct_reads_ignored = 100 * num_ignored / tot_reads if tot_reads > 0 else 0

    # Yell at the user / write out summary stats:
    with open(summary_stats, 'w') as f:
        message = f"Total Reads Processed:\t{tot_reads:d}"
        f.write(f"{message}\n")
        logger.info(message)

        message = f"# Reads Passing Model Filter:\t{num_passed:d}\t{pct_reads_passing:02.4f}"
        f.write(f"{message}\n")
        logger.info(message)

        message = f"# Reads Failing Model Filter:\t{num_failed:d}\t{pct_reads_failing:02.4f}"
        f.write(f"{message}\n")
        logger.info(message)

        message = f"# Reads Ignored:\t{num_ignored:d}\t{pct_reads_ignored:02.4f}"
        f.write(f"{message}\n")
        logger.info(message)

    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)
