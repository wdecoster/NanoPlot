#! /usr/bin/env python
# wdecoster
'''
The main purpose of this script is to create plots for nanopore data.
Input data can be given as
-a directory of basecalled fast5 files
-compressed, standard or streamed fastq file
-compressed, standard or streamed fastq file with additional information added by albacore
-a bam file
For basecalled fast5 files this script optionally performs fastq extraction
'''


from __future__ import division, print_function
import argparse
import sys
import h5py
import os
import time
import logging
import datetime
import re
import seaborn as sns
import pandas as pd
import numpy as np
from multiprocessing import cpu_count
from scipy import stats
import matplotlib.pyplot as plt
import pysam
import nanoget
import nanoplotter
version="0.7.1"

def main():
	'''
	Organization function
	Calls plotting functions with appropriate data
	'''
	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
	stamp = initlogs(time.time())
	datadf, lengthprefix, logBool, stamp = getInput(stamp)
	nanoplotter.scatter(
		x=datadf["readlengths"],
		y=datadf["quals"],
		names=['Read lengths', 'Average read quality'],
		path=os.path.join(args.outdir, args.prefix + lengthprefix + "LengthvsQualityScatterPlot"),
		log=logBool)
	stamp = timeStamp(stamp, "Creating LengthvsQual plot")

	if args.fastq or args.fast5 or args.fastq_albacore:
		nanoplotter.lengthPlots(
			array=datadf["readlengths"],
			name="Read length",
			path=os.path.join(args.outdir, args.prefix + lengthprefix),
			log=logBool)
		stamp = timeStamp(stamp, "Creating length plots")
	if args.fast5 or args.fastq_albacore:
		nanoplotter.spatialHeatmap(
			array=datadf["channelIDs"],
			title="Number of reads gelength outliers for plottingnerated per channel",
			path=os.path.join(args.outdir, args.prefix + "ActivityMap_"),
			filename="ReadsPerChannel",
			colour="Greens")
		stamp = timeStamp(stamp, "Creating spatialheatmap for succesfull basecalls")
		nanoplotter.timePlots(
			df=datadf,
			path=os.path.join(args.outdir, args.prefix))
		stamp = timeStamp(stamp, "Creating timeplots")
	if args.bam:
		bamplots(datadf, lengthprefix, logBool, stamp)
	logging.info("Succesfully processed all input.")


def getInput(stamp):
	'''
	Get input and process accordingly. 	Data can be:
	-a directory of basecalled fast5 files
	-a uncompressed, bgzip, bzip2 or gzip compressed fastq file
	-s sorted bam file
	Handle is passed to the proper functions to get DataFrame with metrics
	'''
	if args.fastq:
		datadf = nanoget.processFastq(args.fastq)
	elif args.fast5:
		if args.fqout:
			sys.stdout = open(os.path.join(args.outdir, args.fqout), 'a')
		elif args.dry:
			sys.stdout = open(os.devnull, 'w')
		datadf, channelfails = nanoget.processFast5(args.fast5, min(cpu_count() - 1, args.threads), args.recursive)
	elif args.bam:
		datadf = nanoget.processBam(args.bam, min(cpu_count() - 1, args.threads))
	elif args.fastq_albacore:
		datadf = nanoget.processFastq_albacore(args.fastq_albacore)
	else:
		logging.error("ARGUMENT ERROR: no input presented.")
		sys.exit('''ARGUMENT ERROR: Required argument is either:\n \
					a (compressed) fastq file [--fastq]\n \
					a directory of basecalled fast5 files [--fast5]\n \
					a bam file [--bam].''')
	stamp = timeStamp(stamp, "Gathering data")
	if args.fast5:
		nanoplotter.spatialHeatmap(
			array=channelfails,
			title="Number of reads which failed basecalling per channel",
			path=os.path.join(args.outdir, args.prefix + "ActivityMap_"),
			filename="BasecallFailedPerChannel",
			colour="OrRd")
		stamp = timeStamp(stamp, "Creating spatialheatmap for basecall fails")
	datadf, lengthprefix, logBool = filterData(datadf)
	return (datadf, lengthprefix, logBool, stamp)


def filterData(datadf):
	'''
	Perform filtering on the data based on arguments set on commandline
	- use alighned length or sequenced length (bam mode only)
	- drop outliers
	- drop reads longer than args.maxlength
	- use log10 scaled reads
	- downsample reads to args.downsample
	Return an accurate prefix which is added to plotnames using this filtered data
	'''
	lengthprefix = []
	if args.alength and args.bam:
		datadf["readlengths"] = datadf["aligned_lengths"]
		lengthprefix.append("Aligned_")
		logging.info("Using aligned read lengths for plotting.")
	else:
		datadf["readlengths"] = datadf["lengths"]
		logging.info("Using sequenced read lengths for plotting.")
	if args.drop_outliers:
		datadf=removeLengthOutliers(datadf, "readlengths")
		lengthprefix.append("OutliersRemoved_")
		logging.info("Removing length outliers for plotting.")
	if args.maxlength:
		datadf=datadf[datadf["readlengths"] < args.maxlength]
		lengthprefix.append("MaxLength" + str(args.maxlength) + '_')
		logging.info("Removing reads longer than {}.".format(str(args.maxlength)))
	if args.loglength:
		datadf["readlengths"] = np.log10(datadf["readlengths"])
		lengthprefix.append("Log_")
		logging.info("Using Log10 scaled read lengths.")
		logBool = True
	else:
		logBool = False
	if args.downsample:
		newNum = min(args.downsample, len(datadf.index))
		datadf = datadf.sample(newNum)
		lengthprefix.append("Downsampled_")
		logging.info("Downsampling the dataset from {} to {} reads".format(len(datadf.index), newNum))
	return(datadf, ''.join(lengthprefix), logBool)


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
	parser.add_argument("--maxlength",
						help="Drop reads longer than length specified.",
						type=int)
	parser.add_argument("--drop_outliers",
						help="Drop outlier reads with extreme long length.",
						action="store_true")
	parser.add_argument("--downsample",
						help="Reduce dataset to N reads by random sampling.",
						type=int)
	parser.add_argument("--loglength",
						help="Logarithmic scaling of lengths in plots.",
						action="store_true")
	parser.add_argument("--alength",
						help="Use aligned read lengths rather than sequenced length (bam mode)",
						action="store_true")
	parser.add_argument("--outdir",
						help="Specify directory in which output has to be created.",
						default=".")
	parser.add_argument("--recursive",
						help="Recursively search the directory for fast5 files.",
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
						help="Data presented is in fastq format.")
	target.add_argument("--fastq_albacore",
						help="Data presented is in fastq format generated by albacore with additional information concerning channel and time.")
	target.add_argument("--fast5",
						help="Data presented in a directory of basecalled fast5 files.")
	target.add_argument("--bam",
						help="Data presented as a bam file.")
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


def initlogs(time0):
	try:
		start_time = datetime.datetime.fromtimestamp(time0).strftime('%Y%m%d_%H%M')
		logging.basicConfig(
			format='%(asctime)s %(message)s',
			filename=os.path.join(args.outdir, args.prefix + "Nanoplot_" + start_time + ".log"),
			level=logging.INFO)
	except IOError:
		sys.exit("ERROR: No writing permission to the directory.")
	logging.info('Nanoplot {} started with arguments {}'.format(version, args))
	logging.info("{} cpu's are available".format(cpu_count()))
	logging.info('Versions of key modules are:')
	for module in [np, h5py, sns, pd, pysam]:
		logging.info('{}: {}'.format(module, module.__version__))
	return time0


def removeLengthOutliers(df, columnname):
	'''Calculation function: Remove records with length-outliers'''
	return df[df[columnname] < (np.median(df[columnname]) + 3 * np.std(df[columnname]))]


def bamplots(datadf, lengthprefix, logBool, stamp):
	'''Call plotting functions specific for bam files'''
	nanoplotter.lengthPlots(
		array=datadf["readlengths"],
		name="Read length",
		path=os.path.join(args.outdir, args.prefix + lengthprefix),
		log=logBool)
	stamp = timeStamp(stamp, "Creating length plots")
	nanoplotter.scatter(
		x=datadf["aligned_lengths"],
		y=datadf["lengths"],
		names=["Aligned read lengths", "Sequenced read length"],
		path=os.path.join(args.outdir, args.prefix + "AlignedReadlengthvsSequencedReadLength"))
	stamp = timeStamp(stamp, "Creating AlignedLengthvsLength plot")
	nanoplotter.scatter(
		x=datadf["mapQ"],
		y=datadf["quals"],
		names=["Read mapping quality", "Average basecall quality"],
		path=os.path.join(args.outdir, args.prefix + "MappingQualityvsAverageBaseQuality"))
	stamp = timeStamp(stamp, "Creating MapQvsBaseQ plot")
	nanoplotter.scatter(
		x=datadf["readlengths"],
		y=datadf["mapQ"],
		names=["Read length", "Read mapping quality"],
		path=os.path.join(args.outdir, args.prefix + lengthprefix + "MappingQualityvsReadLength"),
		log=logBool)
	stamp = timeStamp(stamp, "Creating MapQvsBaseQ plot")
	nanoplotter.scatter(
		x=datadf["percentIdentity"],
		y=datadf["quals"],
		names=["Percent identity", "Read quality"],
		path=os.path.join(args.outdir, args.prefix + "PercentIdentityvsAverageBaseQuality"),
		stat=stats.pearsonr)
	stamp = timeStamp(stamp, "Creating PIDvsBaseQ plot")
	nanoplotter.scatter(
		x=datadf["readlengths"],
		y=datadf["percentIdentity"],
		names=["Aligned read length", "Percent identity"],
		path=os.path.join(args.outdir, args.prefix + "PercentIdentityvsAlignedReadLength"),
		stat=stats.pearsonr,
		log=logBool)
	stamp = timeStamp(stamp, "Creating PIDvsLength plot")


if __name__ == "__main__":
	args = getArgs()
	main()
