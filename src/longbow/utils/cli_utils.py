import functools
import logging
import sys
from pathlib import Path

import click
import click_log

from .constants import DEFAULT_MAX_READ_LENGTH, MODEL_DESC_DELIMITER

log = logging.getLogger(__name__)


def create_logger():
    root_log = logging.getLogger()
    click_log.basic_config(root_log)
    root_log.handlers[0].setFormatter(
        logging.Formatter(
            "[%(levelname)5s %(asctime)s %(name)7s] %(message)s",
            "%Y-%m-%d %H:%M:%S",
        )
    )
    log.debug("Configured base logger")


def input_bam(function):
    return click.argument(
        "input-bam",
        default="-" if not sys.stdin.isatty() else None,
        type=click.File("rb"),
    )(function)


def input_pbi(function):
    return click.option(
        "-p",
        "--pbi",
        required=False,
        type=click.Path(path_type=Path),
        help="BAM .pbi index file",
    )(function)


def output_bam(help_message):
    # note: this one returns a new decorator
    return click.option(
        "-o",
        "--output-bam",
        default="-",
        show_default=True,
        type=click.Path(path_type=Path),
        help=help_message,
    )  # no call here


def reject_bam(function):
    return click.option(
        "-x",
        "--reject-bam",
        default="/dev/null",
        show_default=True,
        type=click.Path(path_type=Path),
        help="Filtered bam output (failing reads only).",
    )(function)


def model(function):
    return click.option(
        "-m",
        "--model",
        help="The model to use for annotation. If not specified, it will be autodetected from "
        f"the BAM header. The value should be of the form [array]{MODEL_DESC_DELIMITER}[cdna] "
        "and both values should be valid models, either preconfigured or loaded by providing "
        "--model-path.",
    )(function)


def force_overwrite(function):
    return click.option(
        "-f",
        "--force",
        is_flag=True,
        default=False,
        show_default=True,
        help="Force overwrite of the output files if they exist.",
    )(function)


def decorator_with_kwarg(function):
    def wrapper(func=None, **kwargs):
        if func is None:
            return functools.partial(function, **kwargs)
        else:
            return function(func, **kwargs)

    return wrapper


@decorator_with_kwarg
def min_length(function, *, additional_help=""):
    return click.option(
        "-l",
        "--min-length",
        type=int,
        default=0,
        show_default=True,
        required=False,
        help="Minimum length of a read to process. Reads shorter than this will not"
        f" be annotated.{additional_help}",
    )(function)


@decorator_with_kwarg
def max_length(function, *, additional_help=""):
    return click.option(
        "-L",
        "--max-length",
        type=int,
        default=DEFAULT_MAX_READ_LENGTH,
        show_default=True,
        required=False,
        help="Maximum length of a read to process. Reads longer than this will not"
        f" be annotated.{additional_help}",
    )(function)


@decorator_with_kwarg
def min_rq(function, *, additional_help=""):
    return click.option(
        "--min-rq",
        type=float,
        default=-2,
        show_default=True,
        required=False,
        help="Minimum ccs-determined read quality for a read to be annotated."
        f" CCS read quality range is [-1,1].{additional_help}",
    )(function)


class MutuallyExclusiveOption(click.Option):
    """Class to define mutually exclusive options for Click.

    -----------------------
    Used ala:
    -----------------------
    @command(help="Run the command.")
    @option('--jar-file', cls=MutuallyExclusiveOption,
            help="The jar file the topology lives in.",
            mutually_exclusive=["other_arg"])
    @option('--other-arg',
            cls=MutuallyExclusiveOption,
            help="The jar file the topology lives in.",
            mutually_exclusive=["jar_file"])
    def cli(jar_file, other_arg):
        print "Running cli."
        print "jar-file: {}".format(jar_file)
        print "other-arg: {}".format(other_arg)

    if __name__ == '__main__':
        cli()


    Taken from: https://stackoverflow.com/a/37491504
    """

    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop("mutually_exclusive", []))
        help_msg = kwargs.get("help", "")
        if self.mutually_exclusive:
            ex_str = ", ".join(self.mutually_exclusive)
            kwargs["help"] = help_msg + (
                " NOTE: This argument is mutually exclusive with "
                " arguments: [" + ex_str + "]."
            )
        super(MutuallyExclusiveOption, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise click.UsageError(
                "Illegal usage: `{}` is mutually exclusive with "
                "arguments `{}`.".format(self.name, ", ".join(self.mutually_exclusive))
            )

        return super(MutuallyExclusiveOption, self).handle_parse_result(ctx, opts, args)


def format_obnoxious_warning_message(message):
    """Adds some obnoxious formatting to the given message."""
    header = r"""
#############################################
__        ___    ____  _   _ ___ _   _  ____
\ \      / / \  |  _ \| \ | |_ _| \ | |/ ___|
 \ \ /\ / / _ \ | |_) |  \| || ||  \| | |  _
  \ V  V / ___ \|  _ <| |\  || || |\  | |_| |
   \_/\_/_/   \_\_| \_\_| \_|___|_| \_|\____|

#############################################

"""

    return f"{header}{message}\n\n#############################################"


def get_field_count_and_percent_string(count, total, fformat="2.4f"):
    count_str = f"{count}/{total}"
    pct_str = f"({100.0*zero_safe_div(count, total):{fformat}}%)"

    return count_str, pct_str


def zero_safe_div(n, d):
    return 0 if not d else n / d
