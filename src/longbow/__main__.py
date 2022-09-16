import logging

import click
import click_log

import importlib
import sys
import pkgutil

import longbow
import longbow.commands
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


# Dynamically find and import sub-commands (allows for plugins at run-time):
for p in pkgutil.iter_modules(longbow.commands.__path__):
    mod = importlib.import_module(f".{p.name}", longbow.commands.__name__)
    main_entry.add_command(getattr(mod, "main"))


if __name__ == "__main__":
    main_entry()  # pylint: disable=E1120
