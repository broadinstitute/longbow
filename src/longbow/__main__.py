import importlib
import logging
import multiprocessing as mp
import pkgutil
import sys

import click
import click_log

import longbow
import longbow.commands

from .meta import VERSION
from .utils.cli_utils import create_logger

logger = logging.getLogger("longbow")


@click.group(name="longbow")
@click_log.simple_verbosity_option(logger, default="WARNING")
@click.option(
    "-t",
    "--threads",
    type=int,
    default=mp.cpu_count() - 1,
    show_default=True,
    help="number of threads to use (0 for all)",
)
@click.pass_context
def main_entry(ctx, threads):
    create_logger()
    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads
    logger.info(f"Running with {threads} worker subprocess(es)")

    ctx.ensure_object(dict)
    ctx.obj["THREADS"] = threads


@main_entry.command()
def version():
    """Print the version of longbow."""
    click.echo(VERSION)


# Dynamically find and import sub-commands (allows for plugins at run-time):
for p in pkgutil.iter_modules(longbow.commands.__path__):
    mod = importlib.import_module(f".{p.name}", longbow.commands.__name__)
    main_entry.add_command(getattr(mod, "main"))


if __name__ == "__main__":
    main_entry()  # pylint: disable=E1120
