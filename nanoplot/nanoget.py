# wdecoster

from __future__ import division
import sys
import time
import logging
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from multiprocessing import Pool
import dateutil.parser
import pysam
import nanoplot.nanomath as nanomath


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
		sys.stderr.write("Terminating worker threads")
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
		quals.append(nanomath.aveQual(read.query_qualities))
		alignedQuals.append(nanomath.aveQual(read.query_alignment_qualities))
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
	channelIDs = []
	for record in SeqIO.parse(inputfastq, "fastq"):
		lengths.append(len(record))
		quals.append(nanomath.aveQual(record.letter_annotations["phred_quality"]))
		channelIDs.append(int(record.description.split('_')[-3].replace('ch', '')))
	datadf["lengths"] = np.array(lengths)
	datadf["quals"] = np.array(quals)
	datadf["channelIDs"] = np.array(channelIDs)
	logging.info("Collected fastq statistics.")
	return datadf


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
	logging.info("Running in fastq mode, expecting additional information added by albacore.")
	inputfastq = handlecompressedFastq(fastq)
	datadf = pd.DataFrame()
	lengths = []
	quals = []
	channels = []
	time_stamps = []
	for record in SeqIO.parse(inputfastq, "fastq"):
		lengths.append(len(record))
		quals.append(nanomath.aveQual(record.letter_annotations["phred_quality"]))
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
