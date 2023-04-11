import logging
import os
import time

import click
import pysam

from ..meta import VERSION
from ..utils import bam_utils, cli_utils

logger = logging.getLogger(__name__)


@click.command("convert")
@cli_utils.output_bam("bam output")
@click.option(
    "-q",
    "--default-rq",
    type=float,
    default=-1,
    show_default=True,
    required=False,
    help="Default read quality for ingested reads (ranging from [-1,1])",
)
@click.option("-g", "--read-group-id", type=str, help="Read group ID to set")
@click.option(
    "-s", "--sample-name", type=str, help="Sample name to set in the read group"
)
@click.option(
    "-l",
    "--min-length",
    type=int,
    default=0,
    required=False,
    help="Minimum length of reads to process",
)
@click.option(
    "-L",
    "--max-length",
    type=int,
    default=30000,
    required=False,
    help="Maximum length of reads to process",
)
@cli_utils.force_overwrite
@click.argument("input-spec", required=True, type=click.Path(exists=True))
def main(
    output_bam,
    default_rq,
    read_group_id,
    sample_name,
    min_length,
    max_length,
    force,
    input_spec,
):
    """Convert reads from fastq{,.gz} files for use with `annotate`."""

    t_start = time.time()

    # Check to see if the output files exist:
    bam_utils.check_for_preexisting_files(output_bam, exist_ok=force)

    # Get the list of input files
    files = _get_files(input_spec)

    # Create header for the output bam file:
    h = pysam.AlignmentHeader.from_dict(
        bam_utils.create_bam_header_with_program_group(
            logger.name,
            pysam.AlignmentHeader().from_dict(
                {
                    "HD": {"VN": f"{VERSION}"},
                    "RG": [{"ID": read_group_id, "SM": sample_name}],
                }
            ),
        )
    )

    num_reads_seen, num_reads_written = 0, 0
    with pysam.AlignmentFile(output_bam, "wb", header=h) as bam_out:
        for file in files:
            with pysam.FastxFile(file) as fh:
                for entry in fh:
                    num_reads_seen += 1

                    if (
                        len(entry.sequence) >= min_length
                        and len(entry.sequence) <= max_length
                    ):
                        r = pysam.AlignedSegment()
                        r.query_name = entry.name
                        r.query_sequence = entry.sequence
                        r.query_qualities = pysam.qualitystring_to_array(entry.quality)
                        r.flag = 4
                        r.mapping_quality = 255
                        r.set_tag("rq", default_rq)
                        r.set_tag("np", 1)
                        r.set_tag("RG", read_group_id)

                        bam_out.write(r)

                        num_reads_written += 1

    et = time.time()
    logger.info(
        f"Done. Elapsed time: {et - t_start:2.2f}s. "
        f"Reads seen: {num_reads_seen}. Reads written: {num_reads_written}."
    )


def _get_files(user_input):
    inputs = []

    if os.path.isdir(user_input):
        for root, dirs, files in os.walk(user_input):
            for file in files:
                if file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
                    inputs.append(os.path.join(root, file))
    elif os.path.isfile(user_input):
        if user_input.endswith(".fastq.gz") or user_input.endswith(".fq.gz"):
            inputs.append(user_input)

    return inputs
