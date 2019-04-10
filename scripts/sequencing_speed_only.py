from nanoget import get_input
from nanoplotter.timeplots import sequencing_speed_over_time
from argparse import ArgumentParser
from os import path


def main():
    args = get_args()
    df = get_input(source="summary", files=args.summary)[["lengths", "start_time", "duration"]]
    sequencing_speed_over_time(dfs=df,
                               path=path.join(args.outdir, args.prefix),
                               figformat=args.format,
                               title=args.title)


def get_args():
    parser = ArgumentParser(
        description="Creates sequencing speed over time plot from summary file.".upper(),
        add_help=False)
    general = parser.add_argument_group(
        title='General options')
    general.add_argument("-h", "--help",
                         action="help",
                         help="show the help and exit")
    general.add_argument("-o", "--outdir",
                         help="Specify directory in which output has to be created.",
                         default=".")
    general.add_argument("-p", "--prefix",
                         help="Specify an optional prefix to be used for the output files.",
                         default="",
                         type=str)
    general.add_argument("--summary",
                         help="Data is in one or more summary file(s) from albacore or guppy.",
                         nargs='+',
                         metavar="file",
                         required=True)
    visual = parser.add_argument_group(
        title='Options for customizing the plots created')
    visual.add_argument("-f", "--format",
                        help="Specify the output format of the plots.",
                        default="png",
                        type=str,
                        choices=['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps',
                                 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff'])
    visual.add_argument("--title",
                        help="Add a title to all plots, requires quoting if using spaces",
                        type=str,
                        default=None)
    return parser.parse_args()


if __name__ == '__main__':
    main()
