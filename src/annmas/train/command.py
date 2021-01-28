import logging
import click
import click_log
import multiprocessing

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command(name="train")
@click_log.simple_verbosity_option(logger)
def main():
    """Train transition and emission probabilities of model on real data"""
    logger.info("train")
