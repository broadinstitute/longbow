import logging
import click
import click_log

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


@click.command(name="inspect")
@click_log.simple_verbosity_option(logger)
def main():
    """Inspect the classification results on specified reads"""
    logger.info("train")
