#!/complgen/bin/anaconda/bin/python2.7
# wdecoster
'''
The main purpose of this script is to create plots for nanopore data.
Input data can be given as
-a directory of basecalled fast5 files
-a directory of raw fast5 files
-compressed, standard or streamed fastq file
-a bam file
For basecalled fast5 files this script also performs fastq extraction
'''

# TODO
# Apply filtering on fast5 to fastq conversion
# Add trimming while writing
# Add adapter detection
# Create one html/pdf/markdown report containing all plots and stats
# Recognize calibration strand DNA, access identity and accuracy
# Asynchronously create plots up to --threads at a time
# Add recursive option for search for fast5 files
# Add option to remove basecalls from fast5, retaining raw reads
# Add optional manual override of fastq path specifying the group
# Overlay heatmaps http://stackoverflow.com/questions/31707033/change-certain-squares-in-a-seaborn-heatmap
# Annotate heatmaps with information about spot on and waste channel http://stackoverflow.com/questions/33158075/custom-annotation-seaborn-heatmap
# Raw fast5 reader for channel and time
# implement optional interactive plots (plotly)
# Does this work with pypy?
# use pylint module
# Stop with sensible error when bam contains no mapped reads
# try-except to catch errors and log the traceback


from __future__ import division, print_function
import argparse
import sys
import h5py
import os
import time
import logging
import datetime
import warnings
import re
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	import seaborn as sns
import pandas as pd
import numpy as np
from multiprocessing import cpu_count
from scipy import stats
import matplotlib.pyplot as plt
# BAD WORKAROUND FOR AVOIDING IMPORT OF WRONG PYSAM MODULE
for dir in ['/storage/breva/complgen/bin/anaconda/lib/python2.7/site-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg', '/complgen/bin/anaconda/lib/python2.7/site-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg']:
	if dir in sys.path:
		sys.path.remove(dir)
import pysam
import nanoget
import nanoplotter
version=0.5.1

def main():
	'''
	Organization function: Get input and process accordingly.
	Data can be:
	-a directory of raw fast5 files
	-a directory of basecalled fast5 files
	-a uncompressed, bgzip, bzip2 or gzip compressed fastq file
	-s sorted bam file
	Handle is passed to the proper functions to get DataFrame with metrics
	Calls plotting functions with appropriate data
	'''
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
	initlogs()
	time0 = time.time()
	if args.fastq:
		if args.fastq == 'stdin':
			logging.info("Reading from stdin.")
			input = sys.stdin  # Expecting uncompressed stream on stdin
		else:  # Annoying that I have to rely on the file extensions
			if args.fastq.endswith('.gz'):
				import gzip
				input = gzip.open(args.fastq, 'r')
			elif args.fastq.endswith('.bz2'):
				import bz2
				input = bz2.BZ2File(args.fastq, 'r')
			elif args.fastq.endswith(('.fastq', '.fq', '.bgz')):
				input = open(args.fastq, 'r')
			else:
				logging.error("INPUT ERROR: Unrecognized file extension")
				sys.exit('''INPUT ERROR: Unrecognized file extension\n,
							supported formats for --fastq are .gz, .bz2, .bgz, .fastq and .fq''')
		datadf = nanoget.processFastq(input)
	elif args.fast5:
		if args.fqout:
			sys.stdout = open(os.path.join(args.outdir, args.fqout), 'a')
		elif args.dry:
			sys.stdout = open(os.devnull, 'w')
		datadf, channelfails = nanoget.processFast5(args.fast5, args.threads, args.recursive)
	elif args.bam:
		datadf = nanoget.processBam(args.bam, args.threads)
	elif args.raw:
		logging.info("Running in raw mode.")
		sys.exit("Not implemented yet!")
	else:
		logging.error("ARGUMENT ERROR: no input presented.")
		sys.exit('''ARGUMENT ERROR: Required argument is either:\n \
					a (compressed) fastq file [--fastq]\n \
					a directory of basecalled fast5 files [--fast5]\n \
					a directory of raw fast5 files [--raw]\n \
					a bam file [--bam].''')
	stamp = timeStamp(time0, "Gathering data")
	if args.downsample:
		prevNum = len(datadf.index)
		newNum = min(10000, len(datadf.index))
		datadf = datadf.sample(newNum)
		logging.info("Downsampled the dataset from {} to {} reads".format(prevNum, newNum))
	#sns.set_style("darkgrid")
	nanoplotter.scatter(
		datadf=datadf,
		var=["lengths", "quals"],
		names=['Read lengths', 'Average read quality'],
		path=os.path.join(args.outdir, args.prefix + "ReadlengthvsQualityScatterPlot"))
	stamp = timeStamp(stamp, "Creating LengthvsQual plot")
	if args.drop_outliers:
		nanoplotter.scatter(
			datadf=removeLengthOutliers(datadf, "lengths"),
			var=["lengths", "quals"],
			names=['Read lengths', 'Average read quality'],
			path=os.path.join(args.outdir, args.prefix + "ReadlengthvsQualityScatterPlot_LengthOutliersRemoved"))
		stamp = timeStamp(stamp, "Creating LengthvsQual plot without length-outliers")
	if args.fastq or args.fast5:
		nanoplotter.lengthPlots(
			array=datadf["lengths"],
			name="Read length",
			path=os.path.join(args.outdir, args.prefix))
		stamp = timeStamp(stamp, "Creating length plots")
		if args.drop_outliers:
			nanoplotter.lengthPlots(
				array=removeLengthOutliers(datadf, "lengths")["lengths"],
				name="Read length",
				path=os.path.join(args.outdir, args.prefix),
				suffix="_LengthOutliersRemoved")
			stamp = timeStamp(stamp, "Creating length plots without length-outliers")
		nanoplotter.spatialHeatmap(
			array=datadf["channelIDs"],
			title="Number of reads generated per channel",
			path=os.path.join(args.outdir, args.prefix + "ActivityMap_"),
			filename="ReadsPerChannel",
			colour="Greens")
		stamp = timeStamp(stamp, "Creating spatialheatmap for succesfull basecalls")
	if args.fast5:
		nanoplotter.spatialHeatmap(
			array=channelfails,
			title="Number of reads which failed basecalling per channel",
			path=os.path.join(args.outdir, args.prefix + "ActivityMap_"),
			filename="BasecallFailedPerChannel",
			colour="OrRd")
		stamp = timeStamp(stamp, "Creating spatialheatmap for basecall fails")
		nanoplotter.timePlots(
			df=datadf,
			path=os.path.join(args.outdir, args.prefix))
		stamp = timeStamp(stamp, "Creating timeplots")
	if args.bam:
		bamplots(datadf, stamp)
	logging.info("Succesfully processed all input.")


def getArgs():
	parser = argparse.ArgumentParser(description="Perform diagnostic plotting, QC analysis and fast5 extraction of Nanopore sequencing data.")
	out = parser.add_mutually_exclusive_group()
	out.add_argument("--dry",
						help="Run on a directory of fast5 files without creating fastq output.",
						action="store_true")
	out.add_argument("--fqout",
						help="File to which output fastq should be written.")
	parser.add_argument("--threads",
						help="Set the allowed number of threads to be used by the script",
						default=4,
						type=int)
	parser.add_argument("--time",
						help="Give timestamps to stderr for optimization purposes",
						action="store_true")
	parser.add_argument("--report",
						help="Summarize all plots in a html report.",
						action="store_true")
	parser.add_argument("--downsample",
						help="Reduce dataset to 10000 reads by random sampling.",
						action="store_true")
	parser.add_argument("--outdir",
						help="Specify directory in which output has to be created.",
						default=".")
	parser.add_argument("--recursive",
						help="Recursively search the directory for fast5 files.",
						action="store_true")
	parser.add_argument("--drop_outliers",
						help="Drop outlier reads with extreme long length.",
						action="store_true")
	parser.add_argument("--version",
						help="Print version and exit.",
						action="store_true")
	parser.add_argument("--prefix",
						help="Specify an optional prefix to be used for the output files.",
						default="",
						type=str)
	target = parser.add_mutually_exclusive_group()
	target.add_argument("--fastq",
						help="Data presented is already in fastq format.")
	target.add_argument("--fast5",
						help="Data presented in a directory of basecalled fast5 files.")
	target.add_argument("--bam",
						help="Data presented as a bam file.")
	target.add_argument("--raw",
						help="Data presented in a directory of raw fast5 files.")
	args = parser.parse_args()
	if args.version:
		print("Nanoplot {}".format(version))
		sys.exit(0)
	return(args)


def timeStamp(start, task):
	now = time.time()
	if args.time:
		sys.stderr.write(task + " took {0:.2f} seconds\n".format(now - start))
	logging.info("Task {0} took {1:.2f} seconds".format(task, now - start))
	return now


def initlogs():
	try:
		start_time = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M')
		logging.basicConfig(
			format='%(asctime)s %(message)s',
			filename=os.path.join(args.outdir, args.prefix + "Nanoplot_" + start_time + ".log"),
			level=logging.INFO)
	except IOError:
		sys.exit("ERROR: No writing permission to the directory.")
	logging.info('Nanoplot started with arguments {}'.format(args))
	logging.info("{} cpu's are available".format(cpu_count()))
	logging.info('Versions of key modules are:')
	for module in [np, h5py, sns, pd, pysam]:
		logging.info('{}: {}'.format(module, module.__version__))


def removeLengthOutliers(df, columnname):
	'''Calculation function: Remove records with length-outliers'''
	return(df[df[columnname] < (np.median(df[columnname]) + 3 * np.std(df[columnname]))])


def removeLowMapQ(df):
	'''Calculation function: Remove records with mapping quality < 10'''
	return(df[df["mapQ"] > 10])


def bamplots(datadf, stamp):
	'''Call plotting functions specific for bam files'''
	nanoplotter.lengthPlots(
		array=datadf["lengths"],
		name="Sequenced read length",
		path=os.path.join(args.outdir, args.prefix))
	nanoplotter.lengthPlots(
		array=datadf["aligned_lengths"],
		name="Aligned read length",
		path=os.path.join(args.outdir, args.prefix))
	stamp = timeStamp(stamp, "Creating length plots")
	nanoplotter.scatter(
		datadf=datadf,
		var=["aligned_lengths", "aligned_quals"],
		names=["Aligned read lengths", "Aligned read quality"],
		path=os.path.join(args.outdir, args.prefix + "AlignedReadlengthvsAlignedQualityScatterPlot"))
	stamp = timeStamp(stamp, "Creating LengthvsQual plot for aligned reads")
	nanoplotter.scatter(
		datadf=datadf,
		var=["aligned_lengths", "lengths"],
		names=["Aligned read lengths", "Sequenced read length"],
		path=os.path.join(args.outdir, args.prefix + "AlignedReadlengthvsSequencedReadLength"))
	stamp = timeStamp(stamp, "Creating AlignedLengthvsLength plot")
	nanoplotter.scatter(
		datadf=datadf,
		var=["mapQ", "quals"],
		names=["Read mapping quality", "Sequenced read quality"],
		path=os.path.join(args.outdir, args.prefix + "MappingQualityvsAverageBaseQuality"))
	stamp = timeStamp(stamp, "Creating MapQvsBaseQ plot")
	nanoplotter.scatter(
		datadf=datadf,
		var=["mapQ", "lengths"],
		names=["Read mapping quality", "Sequenced read length"],
		path=os.path.join(args.outdir, args.prefix + "MappingQualityvsReadLength"))
	stamp = timeStamp(stamp, "Creating MapQvsBaseQ plot")
	nanoplotter.scatter(
		datadf=datadf,
		var=["percentIdentity", "quals"],
		names=["Percent identity", "Sequenced read quality"],
		path=os.path.join(args.outdir, args.prefix + "PercentIdentityvsAverageBaseQuality"),
		stat=stats.pearsonr)
	stamp = timeStamp(stamp, "Creating PIDvsBaseQ plot")
	nanoplotter.scatter(
		datadf=removeLowMapQ(datadf),
		var=["percentIdentity", "quals"],
		names=["Percent identity", "Sequenced read quality"],
		path=os.path.join(args.outdir, args.prefix + "PercentIdentityvsAverageBaseQuality_lowMappingQualityRemoved"),
		stat=stats.pearsonr)
	stamp = timeStamp(stamp, "Creating PIDvsBaseQ plot with low mapping quality reads removed")
	if args.drop_outliers:
		nanoplotter.lengthPlots(
			array=removeLengthOutliers(datadf, "lengths")["lengths"],
			name="Sequenced read length",
			path=os.path.join(args.outdir, args.prefix),
			suffix="_LengthOutliersRemoved")
		nanoplotter.lengthPlots(
			array=removeLengthOutliers(datadf, "lengths")["aligned_lengths"],
			name="Aligned read length",
			path=os.path.join(args.outdir, args.prefix),
			suffix="_LengthOutliersRemoved")
		stamp = timeStamp(stamp, "Creating length plots without length-outliers")
		nanoplotter.scatter(
			datadf=removeLengthOutliers(datadf, "aligned_lengths"),
			var=["aligned_lengths", "aligned_quals"],
			names=["Aligned read lengths", "Aligned read quality"],
			path=os.path.join(args.outdir, args.prefix + "AlignedReadlengthvsAlignedQualityScatterPlot_LengthOutliersRemoved"))
		stamp = timeStamp(stamp, "Creating LengthvsQual plot for aligned reads with length-outliers removed")
		nanoplotter.scatter(
			datadf=removeLengthOutliers(datadf, "lengths"),
			var=["aligned_lengths", "lengths"],
			names=["Aligned read lengths", "Sequenced read length"],
			path=os.path.join(args.outdir, args.prefix + "AlignedReadlengthvsSequencedReadLength_LengthOutliersRemoved"))
		stamp = timeStamp(stamp, "Creating AlignedLengthvsLength plot with length-outliers removed")

if __name__ == "__main__":
	args = getArgs()
	main()
