import importlib
import logging
import multiprocessing as mp
import pkgutil
import sys
from pathlib import Path

import click
import click_log

import longbow
import longbow.commands

from .meta import VERSION
from .utils.cli_utils import create_logger
from .utils.model_utils import load_models

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
@click.option(
    "--model-path",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Path to the model to use for annotation. If a single file, this should "
    "be a JSON file containing both array and cdna structures. If a directory, "
    "the contents of the directory will be loaded and any valid combination of "
    "array and cdna can be specified with --model-name.",
)
@click.pass_context
def main_entry(ctx, threads, model_path=None):
    create_logger()
    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    with importlib.resources.files("longbow.models") as model_dir:
        logger.debug(f"Loading models from {model_dir}")
        load_models(model_dir)

    if model_path is not None:
        logger.debug(f"Loading models from {model_path}")
        load_models(model_path)

    threads = mp.cpu_count() if threads <= 0 or threads > mp.cpu_count() else threads

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
