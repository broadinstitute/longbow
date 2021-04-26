import logging

import click
import click_log

import sys

from .inspect import command as inspect
from .segment import command as segment
from .annotate import command as annotate
from .train import command as train
from .scsplit import command as scsplit
from .discriminate import command as discriminate
from .filter import command as filter

from .meta import VERSION


logger = logging.getLogger("version")
click_log.basic_config(logger)
logger.handlers[0].formatter = logging.Formatter(
    "[%(levelname)s %(asctime)s %(name)8s] %(message)s", "%Y-%m-%d %H:%M:%S"
)


@click.group(name="longbow")
def main_entry():
    logger.info("Invoked via: longbow %s", " ".join(sys.argv))


@main_entry.command()
@click_log.simple_verbosity_option(logger)
def version():
    """Print the version of longbow."""
    click.echo(VERSION)


# Update with new sub-commands:
main_entry.add_command(annotate.main)
main_entry.add_command(segment.main)
main_entry.add_command(train.main)
main_entry.add_command(inspect.main)
main_entry.add_command(scsplit.main)
main_entry.add_command(discriminate.main)
main_entry.add_command(filter.main)


if __name__ == "__main__":
    main_entry()  # pylint: disable=E1120
