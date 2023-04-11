import logging
import sys

import click
import click_log

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
        type=click.Path(exists=False),
        help="BAM .pbi index file",
    )(function)


def output_bam(help_message):
    def decorator(function):
        return click.option(
            "-o",
            "--output-bam",
            default="-",
            show_default=True,
            type=click.Path(exists=False),
            help=help_message,
        )(function)

    return decorator


def reject_bam(function):
    return click.option(
        "-x",
        "--reject-bam",
        default="/dev/null",
        show_default=True,
        type=click.Path(exists=False),
        help="Filtered bam output (failing reads only).",
    )(function)


def model(function):
    return click.option(
        "-m",
        "--model",
        help="The model to use for annotation.  If not specified, it will be autodetected from "
        "the BAM header.  If the given value is a pre-configured model name, then that "
        "model will be used.  Otherwise, the given value will be treated as a file name "
        "and Longbow will attempt to read in the file and create a LibraryModel from it.  "
        "Longbow will assume the contents are the configuration of a LibraryModel as per "
        "LibraryModel.to_json().",
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
        help = kwargs.get("help", "")
        if self.mutually_exclusive:
            ex_str = ", ".join(self.mutually_exclusive)
            kwargs["help"] = help + (
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
