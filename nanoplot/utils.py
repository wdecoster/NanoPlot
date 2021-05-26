import sys
import os
from datetime import datetime as dt
from time import time
import logging
from nanoplot.version import __version__
from argparse import HelpFormatter, Action, ArgumentParser
import textwrap as _textwrap
import pandas as pd


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


class Action_Print_Colormaps(Action):
    def __init__(self, option_strings, dest="==SUPPRESS==", default="==SUPPRESS==", help=None):
        super(Action_Print_Colormaps, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        list_colormaps()


def get_args():
    epilog = """EXAMPLES:
    NanoPlot --summary sequencing_summary.txt --loglength -o summary-plots-log-transformed
    NanoPlot -t 2 --fastq reads1.fastq.gz reads2.fastq.gz --maxlength 40000 --plots hex dot
    NanoPlot --color yellow --bam alignment1.bam alignment2.bam alignment3.bam --downsample 10000
    """
    parser = ArgumentParser(
        description="Creates various plots for long read sequencing data.".upper(),
        epilog=epilog,
        formatter_class=custom_formatter,
        add_help=False)
    general = parser.add_argument_group(
        title='General options')
    general.add_argument("-h", "--help",
                         action="help",
                         help="show the help and exit")
    general.add_argument("-v", "--version",
                         help="Print version and exit.",
                         action="version",
                         version='NanoPlot {}'.format(__version__))
    general.add_argument("-t", "--threads",
                         help="Set the allowed number of threads to be used by the script",
                         default=4,
                         type=int)
    general.add_argument("--verbose",
                         help="Write log messages also to terminal.",
                         action="store_true")
    general.add_argument("--store",
                         help="Store the extracted data in a pickle file for future plotting.",
                         action="store_true")
    general.add_argument("--raw",
                         help="Store the extracted data in tab separated file.",
                         action="store_true")
    general.add_argument("--huge",
                         help="Input data is one very large file.",
                         action="store_true")
    general.add_argument("-o", "--outdir",
                         help="Specify directory in which output has to be created.",
                         default=".")
    general.add_argument("-p", "--prefix",
                         help="Specify an optional prefix to be used for the output files.",
                         default="",
                         type=str)
    general.add_argument("--tsv_stats",
                         help="Output the stats file as a properly formatted TSV.",
                         action='store_true')
    general.add_argument("--info_in_report",
                         help="Add NanoPlot run info in the report.",
                         action='store_true')
    filtering = parser.add_argument_group(
        title='Options for filtering or transforming input prior to plotting')
    filtering.add_argument("--maxlength",
                           help="Hide reads longer than length specified.",
                           type=int,
                           metavar='N')
    filtering.add_argument("--minlength",
                           help="Hide reads shorter than length specified.",
                           type=int,
                           metavar='N')
    filtering.add_argument("--drop_outliers",
                           help="Drop outlier reads with extreme long length.",
                           action="store_true")
    filtering.add_argument("--downsample",
                           help="Reduce dataset to N reads by random sampling.",
                           type=int,
                           metavar='N')
    filtering.add_argument("--loglength",
                           help="Additionally show logarithmic scaling of lengths in plots.",
                           action="store_true")
    filtering.add_argument("--percentqual",
                           help="Use qualities as theoretical percent identities.",
                           action="store_true")
    filtering.add_argument("--alength",
                           help="Use aligned read lengths rather than sequenced length (bam mode)",
                           action="store_true")
    filtering.add_argument("--minqual",
                           help="Drop reads with an average quality lower than specified.",
                           type=int,
                           metavar='N')
    filtering.add_argument("--runtime_until",
                           help="Only take the N first hours of a run",
                           type=int,
                           metavar='N')
    filtering.add_argument("--readtype",
                           help="Which read type to extract information about from summary. \
                                 Options are 1D, 2D, 1D2",
                           default="1D",
                           choices=['1D', '2D', '1D2'])
    filtering.add_argument("--barcoded",
                           help="Use if you want to split the summary file by barcode",
                           action="store_true")
    filtering.add_argument("--no_supplementary",
                           help="Use if you want to remove supplementary alignments",
                           action="store_true",
                           default=False)
    visual = parser.add_argument_group(
        title='Options for customizing the plots created')
    visual.add_argument("-c", "--color",
                        help="Specify a valid matplotlib color for the plots",
                        default="#4CB391")
    visual.add_argument("-cm", "--colormap",
                        help="Specify a valid matplotlib colormap for the heatmap",
                        default="Greens")
    visual.add_argument("-f", "--format",
                        help="Specify the output format of the plots, which are in addition to the html files",
                        default="png",
                        type=str,
                        choices=['png','jpg','jpeg','webp','svg','pdf','eps','json'])
    visual.add_argument("--plots",
                        help="Specify which bivariate plots have to be made.",
                        default=['kde', 'dot'],
                        type=str,
                        nargs='*',
                        choices=['kde', 'hex', 'dot'])
    visual.add_argument("--legacy", help="Specify which bivariate plots have to be made (legacy mode).",
                        type=str,
                        nargs='*',
                        choices=['kde', 'dot', 'hex'])
    visual.add_argument("--listcolors",
                        help="List the colors which are available for plotting and exit.",
                        action=Action_Print_Colors,
                        default=False)
    visual.add_argument("--listcolormaps",
                        help="List the colors which are available for plotting and exit.",
                        action=Action_Print_Colormaps,
                        default=False)
    visual.add_argument("--no-N50",
                        help="Hide the N50 mark in the read length histogram",
                        action="store_true")
    visual.add_argument("--N50",
                        help="Show the N50 mark in the read length histogram",
                        action="store_true",
                        default=False)
    visual.add_argument("--title",
                        help="Add a title to all plots, requires quoting if using spaces",
                        type=str,
                        default=None)
    visual.add_argument("--font_scale",
                        help="Scale the font of the plots by a factor",
                        type=float,
                        default=1)
    visual.add_argument("--dpi",
                        help="Set the dpi for saving images",
                        type=int,
                        default=100)
    visual.add_argument("--hide_stats",
                        help="Not adding Pearson R stats in some bivariate plots",
                        action="store_true",
                        default=False)
    target = parser.add_argument_group(
        title="Input data sources, one of these is required.")
    mtarget = target.add_mutually_exclusive_group(
        required=True)
    mtarget.add_argument("--fastq",
                         help="Data is in one or more default fastq file(s).",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--fasta",
                         help="Data is in one or more fasta file(s).",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--fastq_rich",
                         help="Data is in one or more fastq file(s) generated by albacore, \
                               MinKNOW or guppy with additional information \
                               concerning channel and time.",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--fastq_minimal",
                         help="Data is in one or more fastq file(s) generated by albacore, \
                               MinKNOW or guppy with additional information concerning channel \
                               and time. Is extracted swiftly without elaborate checks.",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--summary",
                         help="Data is in one or more summary file(s) generated by albacore \
                               or guppy.",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--bam",
                         help="Data is in one or more sorted bam file(s).",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--ubam",
                         help="Data is in one or more unmapped bam file(s).",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--cram",
                         help="Data is in one or more sorted cram file(s).",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--pickle",
                         help="Data is a pickle file stored earlier.",
                         metavar="pickle")
    mtarget.add_argument("--feather",
                         help="Data is in one or more feather file(s).",
                         nargs='+',
                         metavar="file")
    args = parser.parse_args()
    if args.listcolors:
        list_colors()
    if args.listcolormaps:
        list_colormaps()
    if args.no_N50:
        sys.stderr.write('DeprecationWarning: --no-N50 is currently the default setting.\n')
        sys.stderr.write('The argument is thus unnecessary but kept for backwards compatibility.')
    if args.barcoded and not args.summary:
        sys.exit('ARGUMENT ERROR: --barcoded only works with data provided as --summary!')
    settings = vars(args)
    settings["path"] = os.path.join(args.outdir, args.prefix)
    return settings, args


def custom_formatter(prog):
    return CustomHelpFormatter(prog)


def list_colors():
    parent_directory = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    colours = open(os.path.join(parent_directory, "extra/color_options_hex.txt"))
    col_hex = {}

    for line in colours:
        key, value = line.split(",")
        col_hex[key] = value.strip()
    print("Valid colors: {}".format("\n".join([c.strip() for c in list(col_hex.keys())])))
    sys.exit(0)


def list_colormaps():
    print('Valid colormaps:\nGreys\nYlGnBu\nGreens\nYlOrRd\nBluered\nRdBu\nReds\nBlues\nPicnic\n'
          'Rainbow\nPortland\nJet\nHot\nBlackbody\nEarth\nElectric\nViridis\nCividis')
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


def subsample_datasets(df, minimal=10000):
    if 'dataset' in df:
        list_df = []

        for d in df["dataset"].unique():
            dataset = df.loc[df['dataset'] == d]

            if len(dataset.index) < minimal:
                list_df.append(dataset)

            else:
                list_df.append(dataset.sample(minimal))

        subsampled_df = pd.concat(list_df, ignore_index=True)

    else:
        if len(df.index) < minimal:
            subsampled_df = df

        else:
            subsampled_df = df.sample(minimal)

    return subsampled_df
