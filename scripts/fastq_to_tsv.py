#! /usr/bin/env python
# wdecoster


import logging
from argparse import ArgumentParser
import sys
import pandas as pd
from Bio import SeqIO
import os
from math import log


def main():
    args = get_args()
    pd.concat(
        [process_fastq_plain(f) for f in args.fastq], ignore_index=True
    ).to_csv(args.output, sep="\t", index=False, compression="gzip")


def handle_compressed_input(inputfq, file_type="fastq"):
    """Return handles from compressed files according to extension.

    Check for which fastq input is presented and open a handle accordingly
    Can read from compressed files (gz, bz2, bgz) or uncompressed
    Relies on file extensions to recognize compression
    """
    if not os.path.isfile(inputfq):
        sys.exit(f"File provided doesn't exist or the path is incorrect: {inputfq}")
    if inputfq.endswith((".gz", "bgz")):
        import gzip

        logging.info("Nanoget: Decompressing gzipped {} {}".format(file_type, inputfq))
        return gzip.open(inputfq, "rt")
    elif inputfq.endswith(".bz2"):
        import bz2

        logging.info("Nanoget: Decompressing bz2 compressed {} {}".format(file_type, inputfq))
        return bz2.open(inputfq, "rt")
    elif inputfq.endswith((".fastq", ".fq", "fasta", ".fa", ".fas")):
        return open(inputfq, "r")
    else:
        logging.error("INPUT ERROR: Unrecognized file extension {}".format(inputfq))
        sys.exit(
            "INPUT ERROR:\nUnrecognized file extension in {}\n"
            "Supported are gz, bz2, bgz, fastq, fq, fasta, fa and fas".format(inputfq)
        )

def process_fastq_plain(fastq):
    """Combine metrics extracted from a fastq file."""
    logging.info("Nanoget: Starting to collect statistics from plain fastq file.")
    inputfastq = handle_compressed_input(fastq)
    return pd.DataFrame(
            data=[res for res in extract_from_fastq(inputfastq) if res],
            columns=["id", "quals", "lengths"],
        ).dropna()


def extract_from_fastq(fq):
    """Extract metrics from a fastq file.

    Return average quality and read length
    """
    for rec in SeqIO.parse(fq, "fastq"):
        yield rec.id, ave_qual(rec.letter_annotations["phred_quality"]), len(rec)


def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10 ** (q / -10) for q in range(n + 1)]


def ave_qual(quals, qround=False, tab=errs_tab(128)):
    """Calculate average basecall quality of a read.

    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        mq = -10 * log(sum([tab[q] for q in quals]) / len(quals), 10)
        if qround:
            return round(mq)
        else:
            return mq
    else:
        return None


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--fastq", help="Data is in one or more default fastq file(s).", nargs="+", metavar="file")
    parser.add_argument("-o", "--output", help="Output raw data to tsv.gz file.", default="fastq_data.tsv.gz")
    return parser.parse_args()


if __name__ == "__main__":
    main()
