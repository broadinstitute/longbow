import logging
import time
import math
import sys

import click
import click_log
import tqdm

import multiprocessing as mp
import concurrent.futures

import pysam

import longbow.utils.constants
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel

from ..utils import bam_utils

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("train")
click_log.basic_config(logger)


@click.command(name="train")
@click_log.simple_verbosity_option(logger)
@click.option(
    "-n",
    "--num-training-samples",
    type=int,
    default=10,
    show_default=True,
    help="number of training samples to use",
)
@click.option(
    "-i",
    "--max-training-iterations",
    type=int,
    default=5,
    show_default=True,
    help="number of training iterations to use",
)
@click.option(
    "-t",
    "--threads",
    type=int,
    default=mp.cpu_count() - 1,
    show_default=True,
    help="number of threads to use (0 for all)",
)
@click.option(
    "-o",
    "--output-yaml",
    required=True,
    type=click.Path(exists=False),
    help="trained model",
)
@click.option(
    "-m",
    "--model",
    default=longbow.utils.constants.DEFAULT_MODEL,
    show_default=True,
    help="The model to use for annotation.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
)
@click.argument("training-bam", type=click.Path(exists=True))
def main(
    num_training_samples, max_training_iterations, threads, output_yaml, model, training_bam
):
    """Train transition and emission probabilities on real data."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    # Get our model:
    if LibraryModel.has_prebuilt_model(model):
        m = LibraryModel.build_pre_configured_model(model)
    else:
        logger.info(f"Loading model from json file: %s", model)
        m = LibraryModel.from_json_file(model)
    logger.info(f"Using %s: %s", model, m.description)

    training_seqs = load_training_seqs(m, num_training_samples, threads, training_bam)

    logger.info("Loaded %d training sequences", len(training_seqs))

    logger.info("Starting training...", len(training_seqs))
    improvement, history = m.fit(
        sequences=training_seqs,
        max_iterations=max_training_iterations,
        stop_threshold=1e-1,
        return_history=True,
        verbose=True,
        n_jobs=threads,
    )

    with open(output_yaml, "w") as model_file:
        print(improvement.to_yaml(), file=model_file)

    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)


def load_training_seqs(m, num_training_samples, threads, training_bam):
    training_seqs = []
    raw_reads = []
    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        training_bam, "rb", check_sq=False, require_index=False
    ) as bam_file:
        for r in bam_file:
            raw_reads.append(r)

            if len(raw_reads) > num_training_samples:
                break
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=threads
    ) as executor, tqdm.tqdm(
        desc="Progress", unit=" reads", colour="green", file=sys.stdout
    ) as pbar:

        future_to_segmented_read = {
            executor.submit(select_read, r, m): r for r in raw_reads
        }

        for future in concurrent.futures.as_completed(future_to_segmented_read):
            read = future_to_segmented_read[future]
            try:
                logp, seq = future.result()
                training_seqs.append(list(seq))
            except Exception as ex:
                logger.error("%r generated an exception: %s", read, ex)

            pbar.update(1)
    return training_seqs


def select_read(read, model):
    flogp = -math.inf
    fseq = None

    # Use the untrained model to determine if we should add this training
    # example in the forward or reverse-complement orientation.
    for seq in [read.query_sequence, bam_utils.reverse_complement(read.query_sequence)]:
        logp, ppath = model.annotate(seq, smooth_islands=True)

        if logp > flogp:
            flogp = logp
            fseq = seq

    return flogp, fseq
