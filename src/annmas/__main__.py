import logging
import click
import click_log

# porcelain
from .segment import command as segment
from .train import command as train
from .inspect import command as inspect

# Version number is automatically set via bumpversion.
# DO NOT MODIFY:
__version__ = "0.0.1"

logger = logging.getLogger(__name__)
click_log.basic_config(logger)
logger.handlers[0].formatter = logging.Formatter("[%(levelname)s %(asctime)s] %(message)s", "%Y-%m-%d %H:%M:%S")


@click.group(name="annmas")
def main_entry():
    pass


@main_entry.command()
@click_log.simple_verbosity_option(logger)
def version():
    """Print the version of annmas"""
    logger.info(f"annmas: {__version__}")


# Update with new sub-commands:
main_entry.add_command(segment.main)
main_entry.add_command(train.main)
main_entry.add_command(inspect.main)


if __name__ == "__main__":
    main_entry()  # pylint: disable=E1120
