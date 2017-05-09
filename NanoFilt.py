#wdecoster

from Bio import SeqIO
import argparse
import sys
from nanoget import aveQualFastq

def getArgs():
    parser = argparse.ArgumentParser(description="Perform quality and or length filtering of Nanopore fastq data on stdin.")
    parser.add_argument("-q", "--quality", help="Filter on a minimum average read quality score", default=0, type=int)
    parser.add_argument("-l", "--length", help="Filter on a minimum read length", default=0, type=int)
    return parser.parse_args()

def filterstream(fq):
    for record in SeqIO.parse(fq, "fastq"):
        if aveQualFastq(record.letter_annotations["phred_quality"]) > args.quality and len(record) > args.length:
            print(record.format("fastq"), end="")

if __name__ == "__main__":
	args = getArgs()
	filterstream(sys.stdin)
