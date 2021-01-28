import logging
import click
import click_log
import sys

import multiprocessing as mp
from multiprocessing import Process, Manager

import concurrent.futures

from collections import namedtuple

from ..utils.model import *


# Named tuple to store alignment information:
class SegmentInfo(namedtuple(
    "SegmentInfo", ["name", "start", "end"]
)):
    def __str__(self):
        return f"SegmentInfo({self.to_tag()})"

    def to_tag(self):
        return f"{self.name}:{self.start}-{self.end}"


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
    logger.info(f"Running with {threads} thread(s)")

    m = build_default_model()
    if model is not None:
        m.from_yaml(model)
        logger.info(f"Using pretrained annotation model {model}")
    else:
        logger.info(f"Using default annotation model")

    num_reads_annotated = 0
    num_sections = 0

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_segmented_read = {executor.submit(segment_read, r, m): r for r in bam_file}

            # Get our header from the input bam file:
            out_bam_header_dict = bam_file.header.to_dict()

            # Add our program group to it:
            pg_dict = {'ID': 'annmas-segment-0.0.1', 'PN': 'annmas', 'VN': '0.0.1',
                       'DS': 'Apply annotation and segmentation model to BAM file.', 'CL': " ".join(sys.argv)}
            if "PG" in out_bam_header_dict:
                out_bam_header_dict["PG"].append(pg_dict)
            else:
                out_bam_header_dict["PG"] = [pg_dict]

            with pysam.AlignmentFile(output_bam, 'wb', header=out_bam_header_dict) as out_bam_file:

                for future in concurrent.futures.as_completed(future_to_segmented_read):
                    read = future_to_segmented_read[future]
                    try:
                        path, logp = future.result()
                    except Exception as ex:
                        logger.error('%r generated an exception: %s', read, ex)
                    else:
                        # Condense the output annotations so we can write them out with indices:
                        segments = collapse_annotations(path)

                        # Obligatory log message:
                        logger.debug("Path for read %s (%2.2f): %s", read.query_name, logp, segments)

                        # Set our tag and write out the read:
                        read.set_tag("SG", "|".join([s.to_tag() for s in segments]))
                        out_bam_file.write(read)

                        # Increment our counters:
                        num_reads_annotated += 1
                        num_sections += len(segments)

    logger.info("annmas: segment finished")
    logger.info(f"annmas: segmented {num_reads_annotated} reads into {num_sections} sections")


def collapse_annotations(path):
    last = ""
    start = 0
    segments = []
    for i, seg in enumerate(path):
        if seg != last:
            if i != 0:
                segments.append(SegmentInfo(last, start, i - 1))
            last = seg
            start = i
    # Don't forget the last one:
    segments.append(SegmentInfo(last, start, i - 1))
    return segments


def segment_read(read, um):
    flogp = -math.inf
    for seq in [read.query_sequence, reverse_complement(read.query_sequence)]:
        logp, ppath = annotate(um, seq)

        if logp > flogp:
            flogp = logp
            fppath = ppath

    return fppath, flogp
