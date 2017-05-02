# NanoPlot
Plotting scripts for nanopore sequencing data

This script performs data extraction from Nanopore sequencing data in the following formats:
-directories of fast5 files
-bam files
-fastq files (can be bgzip, bzip2 or gzip compressed)

Fastq can be extracted while getting data from fast5 files, streaming to stdout or to a file.

Various plots are created automatically, some optionally.

## USAGE:
usage: NanoPlot.py [-h] [--dry | --fqout FQOUT] [--threads THREADS] [--time]
                   [--report] [--downsample] [--outdir OUTDIR] [--recursive]
                   [--drop_outliers] [--prefix PREFIX]
                   [--fastq FASTQ | --fast5 FAST5 | --bam BAM | --raw RAW]

Perform diagnostic plotting, QC analysis and fast5 extraction of Nanopore
sequencing data.

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
