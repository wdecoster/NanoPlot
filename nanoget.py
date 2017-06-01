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
	if not samfile.header['HD']['SO'] == 'coordinate':
		logging.info("Bam file not sorted by coordinate!.")
		sys.exit("Please use a bam file sorted by coordinate.")
	NumberOfmappedReads = samfile.mapped
	NumberOfunmappedReads = samfile.unmapped
	logging.info("Bam file contains {} mapped and {} unmapped reads.".format(NumberOfmappedReads, NumberOfunmappedReads))
	if NumberOfmappedReads == 0:
		sys.exit("FATAL: not a single read was mapped in the bam file.")
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
	datadf["percentIdentity"] = np.array([x for y in [elem[5] for elem in output] for x in y])
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
	pID = []
	for read in samfile.fetch(reference=chromosome, multiple_iterators=True):
		lengths.append(read.query_length)
		alignedLengths.append(read.query_alignment_length)
		quals.append(aveQualBam(read.query_qualities))
		alignedQuals.append(aveQualBam(read.query_alignment_qualities))
		mapQ.append(read.mapping_quality)
		try:
			pID.append((1- read.get_tag("NM")/read.query_alignment_length)*100)
		except KeyError:
			pID.append((1 -
				( parseMD(read.get_tag("MD")) + parseCIGAR(read.cigartuples))
				/read.query_alignment_length)*100)
	return (lengths, alignedLengths, quals, alignedQuals, mapQ, pID)

def parseMD(MDlist):
	return sum([len(item) for item in re.split('[0-9^]', MDlist )])


def parseCIGAR(cigartuples):
	return sum([item[1] for item in cigartuples if item[0] == 1])


def aveQualBam(quals):
	return sum(quals) / len(quals)


def processFast5(directory, threads, recursive):
	'''
	Processing function, calls worker function
	Organize processing basecalled fast5 data, extraction and gathering statistics
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
					print(fq)  # Note that stdout gets redirected with args.dry or args.fqout
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
	a_time_stamps = np.array(times, dtype='datetime64[s]')
	datadf["start_time"] = (a_time_stamps - np.amin(a_time_stamps))
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


def handlecompressedFastq(inputfq):
	'''
	Check for which fastq input is presented and open a handle accordingly
	Can read from stdin, compressed files (gz, bz2, bgz) or uncompressed
	Relies on file extensions to recognize compression
	'''
	if inputfq == 'stdin':
		logging.info("Reading from stdin.")
		return sys.stdin
	else:
		if inputfq.endswith('.gz'):
			import gzip
			logging.info("Decompressing gzipped fastq.")
			return gzip.open(inputfq, 'rt')
		elif inputfq.endswith('.bz2'):
			import bz2
			logging.info("Decompressing bz2 compressed fastq.")
			return bz2.BZ2File(inputfq, 'rt')
		elif inputfq.endswith(('.fastq', '.fq', '.bgz')):
			return open(inputfq, 'r')
		else:
			logging.error("INPUT ERROR: Unrecognized file extension")
			sys.exit('''INPUT ERROR: Unrecognized file extension\n,
						supported formats for --fastq are .gz, .bz2, .bgz, .fastq and .fq''')

def processFastq(fastq):
	'''
	Processing function
	Iterate over a fastq file and extract metrics
	'''
	logging.info("Running in fastq mode.")
	inputfastq = handlecompressedFastq(fastq)
	datadf = pd.DataFrame()
	lengths = []
	quals = []
	for record in SeqIO.parse(inputfastq, "fastq"):
		lengths.append(len(record))
		quals.append(aveQualFastq(record.letter_annotations["phred_quality"]))
	datadf["lengths"] = np.array(lengths)
	datadf["quals"] = np.array(quals)
	logging.info("Collected fastq statistics.")
	return datadf


def aveQualFastq(quals):
	'''	Calculation function: Receive the integer quality scores of a read and return the average quality for that read'''
	return sum(quals) / len(quals)


def processFastq_albacore(fastq):
	'''
	Extract information from fastq files generated by albacore, containing richer information in the header
	containing key-value pairs
	read=<int> [72]
	ch=<int> [159]
	start_time=<timestamp> [2016-07-15T14:23:22Z]  # UTC ISO 8601 ISO 3339 timestamp
	Z indicates UTC time, T is the delimiter between date expression and time expression
	dateutil.parser.parse("2016-07-15T14:23:22Z") # -> datetime.datetime(2016, 7, 15, 14, 23, 22, tzinfo=tzutc())
	'''
	import dateutil.parser
	logging.info("Running in fastq mode, expecting additional information added by albacore.")
	inputfastq = handlecompressedFastq(fastq)
	datadf = pd.DataFrame()
	lengths = []
	quals = []
	channels = []
	time_stamps = []
	for record in SeqIO.parse(inputfastq, "fastq"):
		lengths.append(len(record))
		quals.append(aveQualFastq(record.letter_annotations["phred_quality"]))
		for data in record.description.split(' '):  # This can easily be adapted to include more metrics using the same format
			if data.startswith('ch='):
				channels.append(int(data[3:]))
			elif data.startswith('start_time='):
				time_stamps.append(dateutil.parser.parse(data[11:]))
	datadf["lengths"] = np.array(lengths)
	datadf["quals"] = np.array(quals)
	datadf["channelIDs"] = np.array(channels)
	a_time_stamps = np.array(time_stamps, dtype='datetime64[s]')
	datadf["start_time"] = a_time_stamps - np.amin(a_time_stamps)
	logging.info("Collected fastq statistics.")
	return datadf
