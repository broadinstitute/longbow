import logging
import time
import os
import sys
from collections import Counter

import click

import tqdm
import pysam

import longbow.utils.constants
from ..utils import bam_utils
from ..utils import cli_utils
from ..utils.bam_utils import get_segments
from ..utils.constants import FFORMAT
from ..utils.cli_utils import zero_safe_div

logger = logging.getLogger(__name__)


@click.command()
@cli_utils.input_pbi
@cli_utils.output_bam("filtered bam output (passing reads only)")
@click.option(
    "-x",
    "--reject-bam",
    default="/dev/null",
    type=click.Path(exists=False),
    help="Filtered bam output (failing reads only)  [default: /dev/null]",
)
@cli_utils.model
@cli_utils.force_overwrite
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
@cli_utils.input_bam
def main(pbi, output_bam, reject_bam, model, force, stats, summary_stats, ignore_list, input_bam):
    """Filter segmented reads by conformation to expected cDNA design."""

    t_start = time.time()

    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("About to Sift %d reads", read_count)
    else:
        read_count = bam_utils.get_read_count_from_bam_index(input_bam)
        if read_count:
            logger.info("About to Sift %d reads", read_count)

    # Get our model:
    lb_model = bam_utils.load_model(model, input_bam)
    logger.info(f"Using {lb_model.name}: {lb_model.description}")

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
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file:

        # Get our header from the input bam file:
        out_header = pysam.AlignmentHeader.from_dict(
            bam_utils.create_bam_header_with_program_group(logger.name, bam_file.header)
        )

        # Define a hard print interval that will periodically let the user know something is going on:
        if read_count:
            hard_print_interval = int(read_count * 0.1)
        else:
            hard_print_interval = 100000

        # Setup output files:
        with pysam.AlignmentFile(output_bam, "wb", header=out_header) as passing_bam_file, \
                pysam.AlignmentFile(reject_bam, "wb", header=out_header) as failing_bam_file, \
                open(stats, 'w') as stats_file:

            num_passed = 0
            num_failed = 0
            num_ignored = 0

            all_model_states = sorted(lb_model.cdna_model['structure'] + lb_model.key_adapters + ["random"])

            stats_file.write('\t'.join(all_model_states) + '\n')

            for i, read in enumerate(tqdm.tqdm(bam_file, desc="Progress", unit=" read", colour="green", file=sys.stderr,
                                               disable=not sys.stdin.isatty(), total=read_count)):

                if read.query_name in reads_to_ignore:
                    logger.debug(f"Ignoring read: {read.query_name}")
                    num_ignored += 1
                    continue

                # Get our read segments:
                try:
                    seq, segment_ranges, segment_cigars = get_segments(read)
                except KeyError:
                    logger.error(f"Input bam file does not contain longbow segmented reads!  "
                                 f"No {longbow.utils.constants.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                    sys.exit(1)

                is_valid, actual_element_counts = check_validity(lb_model, segment_ranges)

                if is_valid:
                    logger.debug("Read is %s valid: %s",
                                 lb_model.name,
                                 read.query_name)

                    passing_bam_file.write(read)
                    num_passed += 1
                else:
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug("Read is not %s valid: %s",
                                     lb_model.name,
                                     read.query_name)

                    failing_bam_file.write(read)

                    # Write out read info to stats file:
                    stats_file.write('\t'.join(
                        [str(actual_element_counts[s]) for s in all_model_states]
                    ) + '\n')

                    num_failed += 1

                if (i % hard_print_interval) == 0:
                    if read_count:
                        logger.info("Sifted %d/%d (%2.2f%%) reads.", i, read_count, 100*(i/read_count))
                    else:
                        logger.info("Sifted %d reads.", i)

    # Calc some stats:
    tot_reads = num_passed + num_failed + num_ignored

    # Yell at the user / write out summary stats:
    with open(summary_stats, 'w') as f:
        message = f"Total Reads Processed:\t{tot_reads:d}"
        f.write(f"{message}\n")
        logger.info(message)

        message = f"# Reads Passing Model Filter:\t{num_passed:d}\t{100*zero_safe_div(num_passed, tot_reads):{FFORMAT}}"
        f.write(f"{message}\n")
        logger.info(message)

        message = f"# Reads Failing Model Filter:\t{num_failed:d}\t{100*zero_safe_div(num_failed, tot_reads):{FFORMAT}}"
        f.write(f"{message}\n")
        logger.info(message)

        message = f"# Reads Ignored:\t{num_ignored:d}\t{100*zero_safe_div(num_ignored, tot_reads):{FFORMAT}}"
        f.write(f"{message}\n")
        logger.info(message)

    logger.info(f"Done. Elapsed time: %{FFORMAT}s.", time.time() - t_start)


def check_validity(lb_model, segment_ranges):
    expected_elements = lb_model.cdna_model['structure']

    actual_elements = []
    for s in segment_ranges:
        actual_elements.append(s.name)

    valid = True
    if expected_elements[0] in actual_elements:
        i = actual_elements.index(expected_elements[0])

        for j in range(len(expected_elements)):
            valid &= i+j < len(actual_elements) and expected_elements[j] == actual_elements[i+j]

    valid &= 'random' not in actual_elements

    c = Counter(actual_elements)
    for e in expected_elements:
        valid &= c[e] == 1

    return valid, c
