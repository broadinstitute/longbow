import logging
import click
import click_log

import multiprocessing as mp
from multiprocessing import Process, Manager

from ..utils.model import *

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command(name="segment")
@click_log.simple_verbosity_option(logger)
@click.option("-m", "--model", required=False, type=click.Path(exists=True), help="pre-trained model to apply")
@click.option("-t", "--threads", type=int, default=1, show_default=True, help="number of threads to use (0 for all)")
@click.option("-o", "--output-bam", required=True, type=click.Path(exists=False), help="segmented bam output")
@click.argument('input-bam', type=click.Path(exists=True))
def main(model, threads, output_bam, input_bam):
    """Apply annotation and segmentation model to BAM file"""
    logger.info(f"annmas: segment started")

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"running with {threads} thread(s)")

    m = build_default_model()
    if model is not None:
        m.from_yaml(model)
        logger.info(f"using pretrained annotation model {model}")
    else:
        logger.info(f"using default annotation model")

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bf:
        for read in bf:
            print(read.query_name)
            segment_read(read, m)

    logger.info("annmas: segment finished")


def segment_read(read, um):
    flogp = -math.inf
    for seq in [read.query_sequence, reverse_complement(read.query_sequence)]:
        logp, ppath = annotate(um, seq)

        if logp > flogp:
            flogp = logp
            fppath = ppath
