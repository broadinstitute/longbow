import logging

import click
import click_log

import sys
import pkgutil

import longbow
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
# An alternative would be to iterate through the following:
# pkgutil.iter_modules([os.path.dirname((inspect.getmodule(sys.modules[__name__]).__file__))])
for p in [p for p in pkgutil.iter_modules(longbow.__path__) if p.ispkg]:
    exec(f"from .{p.name} import command as {p.name}")
    exec(f"main_entry.add_command({p.name}.main)")


if __name__ == "__main__":
    main_entry()  # pylint: disable=E1120
