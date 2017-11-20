import sys
import os
from datetime import datetime as dt
from time import time
import logging
from nanoplot.version import __version__
from argparse import HelpFormatter, Action
import textwrap as _textwrap


class CustomHelpFormatter(HelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string

    def _fill_text(self, text, width, indent):
        return ''.join(indent + line for line in text.splitlines(keepends=True))

    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, 80)


class Action_Print_Colors(Action):
    def __init__(self, option_strings, dest="==SUPPRESS==", default="==SUPPRESS==", help=None):
        super(Action_Print_Colors, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        list_colors()


def custom_formatter(prog):
    return CustomHelpFormatter(prog)


def list_colors():
    parent_directory = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    colours = open(os.path.join(parent_directory, "extra/color_options.txt")).readlines()
    print("{}".format(", ".join([c.strip() for c in colours])))
    sys.exit(0)


def make_output_dir(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except IOError:
        sys.exit("ERROR: No writing permission to the output directory.")


def init_logs(args, tool="NanoPlot"):
    """Initiate log file and log arguments."""
    start_time = dt.fromtimestamp(time()).strftime('%Y%m%d_%H%M')
    logname = os.path.join(args.outdir, args.prefix + tool + "_" + start_time + ".log")
    handlers = [logging.FileHandler(logname)]
    if args.verbose:
        handlers.append(logging.StreamHandler())
    logging.basicConfig(
        format='%(asctime)s %(message)s',
        handlers=handlers,
        level=logging.INFO)
    logging.info('{} {} started with arguments {}'.format(tool, __version__, args))
    logging.info('Python version is: {}'.format(sys.version.replace('\n', ' ')))
    return logname
