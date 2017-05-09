# NanoPlot
Plotting scripts for nanopore sequencing data

This script performs data extraction from Nanopore sequencing data in the following formats:
- directories of fast5 files  
- bam files  
- fastq files (can be bgzip, bzip2 or gzip compressed)  

Fastq can be extracted while getting data from fast5 files, streaming to stdout or to a file.

Various plots are created automatically, some optionally.

## USAGE:
usage: NanoPlot.py [-h] [--dry | --fqout FQOUT] [--threads THREADS] [--time]
                   [--report] [--downsample] [--outdir OUTDIR] [--recursive]
                   [--drop_outliers] [--prefix PREFIX]
                   [--fastq FASTQ | --fast5 FAST5 | --bam BAM | --raw RAW]


optional arguments:  
  -h, --help         show this help message and exit  
  --dry              Run on a directory of fast5 files without creating fastq output.  
  --fqout FQOUT      File to which output fastq should be written.  
  --threads THREADS  Set the allowed number of threads to be used by the script  
  --time             Give timestamps to stderr for optimization purposes  
  --downsample       Reduce dataset to 10000 reads by random sampling.  
  --outdir OUTDIR    Specify directory in which output has to be created.  
  --recursive        Recursively search the directory for fast5 files.  
  --drop_outliers    Drop outlier reads with extreme long length.  
  --prefix PREFIX    Specify an optional prefix to be used for the output files.  
  --fastq FASTQ      Data presented is already in fastq format.  
  --fast5 FAST5      Data presented in a directory of basecalled fast5 files.  
  --bam BAM          Data presented as a bam file.  

# Companion script NanoFilt.py
Script to perform filtering on quality and/or read length, and optional trimming after passing filters
Reads from stdin, writes to stdout
Intended to be used:
-directly after fastq extraction
-prior to mapping
-in a stream between extraction and mapping

## USAGE:
usage: NanoFilt.py [-h] [-q QUALITY] [-l LENGTH]

optional arguments:  
  -h, --help            show this help message and exit  
  -q QUALITY, --quality QUALITY  Filter on a minimum average read quality score  
  -l LENGTH, --length LENGTH Filter on a minimum read length  
  --headcrop HEADCROP   Trim n nucleotides from start of read  
  --tailcrop TAILCROP   Trim n nucleotides from end of read
