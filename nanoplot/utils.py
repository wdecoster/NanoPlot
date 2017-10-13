import sys
import os
from datetime import datetime as dt
from time import time
import logging
from nanoplot.version import __version__


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
    '''
    Initiate log file
    Log arguments and module versions
    '''
    start_time = dt.fromtimestamp(time()).strftime('%Y%m%d_%H%M')
    logging.basicConfig(
        format='%(asctime)s %(message)s',
        filename=os.path.join(args.outdir, args.prefix + tool + "_" + start_time + ".log"),
        level=logging.INFO)
    logging.info('{} {} started with arguments {}'.format(tool, __version__, args))
    logging.info('Python version is: {}'.format(sys.version.replace('\n', ' ')))
