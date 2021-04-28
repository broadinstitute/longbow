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

# Update with new sub-commands:
# NOTE: Keeping ths around in case we hate my solution above.

# from .inspect import command as inspect
# from .segment import command as segment
# from .annotate import command as annotate
# from .train import command as train
# from .scsplit import command as scsplit
# from .discriminate import command as discriminate
# from .filter import command as filter
# from .extract import command as extract

# main_entry.add_command(annotate.main)
# main_entry.add_command(segment.main)
# main_entry.add_command(train.main)
# main_entry.add_command(inspect.main)
# main_entry.add_command(scsplit.main)
# main_entry.add_command(discriminate.main)
# main_entry.add_command(filter.main)
# main_entry.add_command(extract.main)

if __name__ == "__main__":
    main_entry()  # pylint: disable=E1120
