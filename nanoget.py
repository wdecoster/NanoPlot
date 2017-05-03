#!/complgen/bin/anaconda/bin/python2.7
# wdecoster

from __future__ import division, print_function
import sys
import h5py
import os
import time
import logging
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from multiprocessing import Pool
# BAD WORKAROUND FOR AVOIDING IMPORT OF WRONG PYSAM MODULE
for dir in ['/storage/breva/complgen/bin/anaconda/lib/python2.7/site-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg', '/complgen/bin/anaconda/lib/python2.7/site-packages/RSeQC-2.6.3-py2.7-linux-x86_64.egg']:
	if dir in sys.path:
		sys.path.remove(dir)
import pysam


def processBam(bam, threads):
	'''
	Processing function: calls pool of worker functions
	to extract from a bam file the following metrics:
	-lengths
	-aligned lengths
	-qualities
	-aligned qualities
	-mapping qualities
	-edit distances to the reference genome scaled by read length
	Returned in a pandas DataFrame
	'''
	logging.info("Running in bam mode.")
	samfile = pysam.AlignmentFile(bam, "rb")
	if not samfile.has_index():
		pysam.index(bam)
		samfile = pysam.AlignmentFile(bam, "rb")  # Need to reload the samfile after creating index
		logging.info("No index for bam file could be found, created index.")
	NumberOfmappedReads = samfile.mapped
	NumberOfunmappedReads = samfile.unmapped
	logging.info("Bam file contains {} mapped and {} unmapped reads.".format(NumberOfmappedReads, NumberOfunmappedReads))
	chromosomes = samfile.references
	datadf = pd.DataFrame()
	pool = Pool(processes=threads)
	params = zip([bam]*len(chromosomes), chromosomes)
	try:
		output = [results for results in pool.imap(extractFromBam, params)]
	except KeyboardInterrupt:
		print("Terminating worker threads")
		pool.terminate()
		pool.join()
		sys.exit()
	# Output contains a tuple per worker
	# Each tuple contains lists per metric
	# Unpacked by following nested list comprehensions
	datadf["lengths"] = np.array([x for y in [elem[0] for elem in output] for x in y])
	datadf["aligned_lengths"] = np.array([x for y in [elem[1] for elem in output] for x in y])
	datadf["quals"] = np.array([x for y in [elem[2] for elem in output] for x in y])
	datadf["aligned_quals"] = np.array([x for y in [elem[3] for elem in output] for x in y])
	datadf["mapQ"] = np.array([x for y in [elem[4] for elem in output] for x in y])
	datadf["editDistances"] = np.array([x for y in [elem[5] for elem in output] for x in y])
	assert datadf["lengths"].size == NumberOfmappedReads, "Unexpected difference in length of entries in datadict"
	logging.info("Collected bam statistics.")
	return datadf


def extractFromBam(params):
	'''
	Worker function per chromosome
	loop over a bam file and create tuple with lists containing metrics:
	-lengths
	-aligned lengths
	-qualities
	-aligned qualities
	-mapping qualities
	-edit distances to the reference genome scaled by read length
	'''
	bam, chromosome = params
	samfile = pysam.AlignmentFile(bam, "rb")
	lengths = []
	alignedLengths = []
	quals = []
	alignedQuals = []
	mapQ = []
	editDistances = []
	for read in samfile.fetch(reference=chromosome, multiple_iterators=True):
		lengths.append(read.query_length)
		alignedLengths.append(read.query_alignment_length)
		quals.append(aveQualBam(read.query_qualities))
		alignedQuals.append(aveQualBam(read.query_alignment_qualities))
		mapQ.append(read.mapping_quality)
		try:
			editDistances.append(read.get_tag("NM")/read.query_length)
		except KeyError:
			editDistances.append(
				(sum([len(item) for item in re.split('[0-9^]', read.get_tag("MD"))]) +  # Parse MD string to get mismatches/deletions
				sum([item[1] for item in read.cigartuples if item[0] == 1]))  # Parse cigar to get insertions
				/read.query_alignment_length)
	return (lengths, alignedLengths, quals, alignedQuals, mapQ, editDistances)


def aveQualBam(quals):
	return sum(quals) / len(quals)


def processFast5(directory, threads, recursive):
	'''
	Processing function, calls worker function
	Organize processing fast5 data, extraction and gathering statistics
	'''
	logging.info("Running in fast5 mode.")
	if recursive:
		fast5list = [os.path.join(root,f) for root, __, files in os.walk(directory) for f in files if f.endswith(".fast5")]
	else:
		fast5list = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".fast5")]
	datadf = pd.DataFrame()
	lengths = []
	quals = []
	channels = []
	channelfails = []
	times = []
	if not fast5list:
		logging.error("No fast5 file found in specified directory {}.".format(directory))
		sys.exit("Not a single fast5 file found in specified directory {}".format(directory))
	else:
		logging.info("Found {} fast5 files.".format(len(fast5list)))
		fqpath, timepath = tasteFast5(fast5list[:101])
		logging.info("Set path to fastq data to {}".format(fqpath))
		pool = Pool(processes=threads)
		params = zip(fast5list, [fqpath] * len(fast5list), [timepath] * len(fast5list))
		logging.info("Collecting fast5 statistics, using {} threads".format(threads))
		try:
			for results in pool.imap(extractFromFast5, params):
				if isinstance(results, tuple):
					l, q, cID, fq, t = results
					lengths.append(l)
					quals.append(q)
					channels.append(cID)
					times.append(t)
					print(fq)
				else:
					channelfails.append(results)
		except KeyboardInterrupt:
				print("Terminating worker threads")
				pool.terminate()
				pool.join()
				sys.exit()
	datadf["lengths"] = np.array(lengths)
	datadf["quals"] = np.array(quals)
	datadf["channelIDs"] = np.array(channels)
	datadf["start_time"] = np.array(times)
	logging.info("Collected fast5 statistics.")
	logging.info("No basecall was made in {} fast5 files".format(len(channelfails)))
	return (datadf, np.array(channelfails))


def tasteFast5(fast5files):
	'''
	Helper function
	Will open fast5 files and find the path to the most recent basecall.
	As soon as a path is found the path is returned
	This assumes all files are from the same dataset.
	If no path is found in subset of files given (100), give error and stop
	Got some ideas from https://github.com/rrwick/Fast5-to-Fastq/blob/master/fast5_to_fastq.py
	'''
	for fast5file in fast5files:
		try:
			with h5py.File(fast5file, 'r') as hdf:
				names = get_hdf5_names(hdf)
				latestFQ = sorted([x for x in names if x.upper().endswith('FASTQ')])[-1]
				#latestCall = [entry for entry in hdf["/Analyses"].keys() if entry.startswith('Basecall')][-1]
				logging.info("Established {} as latest basecall in fast5 files".format(latestFQ))
				timepath = ["/Analyses/" + entry + "/BaseCalled_template/Events" for entry in hdf["/Analyses"].keys() if entry.startswith('Basecall') and "start_time" in hdf["/Analyses/" + entry + "/BaseCalled_template/Events"].attrs.keys()][0]
				logging.info("Established {} as path to creation time of fast5 file".format(timepath))
				return(latestFQ, timepath)
		except KeyError:
			logging.info("Encountered a file lacking basecall, trying the next to find the path to the latest basecalls.")
			continue
	else:
		logging.error("No basecall was found in the first 100 reads.")
		sys.exit("ERROR: a basecall could not be found in the first 100 reads.")


def get_hdf5_names(hdf5_file):
    names = []
    hdf5_file.visit(names.append)
    return names


def extractFromFast5(params):
	'''
	Worker function
	Extracts fastq data from fast5 file, based on fqpath
	Could break for older versions of the fast5 format!
	fastq is printed to stdout,
	aveQualFast5() gets data and returns tuple of length, average quality and channel_number
	'''
	fast5file, fqpath, timepath = params
	with h5py.File(fast5file, 'r') as hdf:
		try:
			return(aveQualFast5(
				hdf[fqpath][()].rstrip(),
				hdf["/UniqueGlobalKey/channel_id"].attrs['channel_number'],
				hdf[timepath].attrs['start_time']))
		except KeyError:
			return(int(hdf["/UniqueGlobalKey/channel_id"].attrs['channel_number']))


def aveQualFast5(fastq, channelID, time):
	'''
	Calculation function
	Receive the phred ascii quality scores of a read and return the length and average float quality for that read
	ASCII codes are !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
	which translate to integers 0-93 by getting ord() of ascii minus 33
	the channelID and time are simply here to get packed in the returned tuple
	'''
	quals = fastq.split(b'\n')[3].decode("utf-8")
	return(
		len(quals),
		sum([ord(elem) - 33 for elem in quals]) / len(quals),
		int(channelID),
		fastq.decode("utf-8"),
		int(time))


def processFastq(inputfastq):
	'''
	Processing function
	Iterate over a fastq file and extract metrics
	'''
	logging.info("Running in fastq mode.")
	datadf = pd.DataFrame()
	lengths = []
	quals = []
	channels = []
	for record in SeqIO.parse(inputfastq, "fastq"):
		lengths.append(len(record))
		quals.append(aveQualFastq(record.letter_annotations["phred_quality"]))
		channels.append(int(record.description.split('_')[-3][2:]))
	datadf["lengths"] = np.array(lengths)
	datadf["quals"] = np.array(quals)
	datadf["channelIDs"] = np.array(channels)
	logging.info("Collected fastq statistics.")
	return datadf


def aveQualFastq(quals):
	'''	Calculation function: Receive the integer quality scores of a read and return the average quality for that read'''
	return sum(quals) / len(quals)
