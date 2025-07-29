#!/bin/bash
set -ev

if [ -d "nanotest" ]; then
    echo "nanotest already cloned"
else
    git clone https://github.com/wdecoster/nanotest.git
fi

# Quick info tests (no data processing)
NanoPlot -h
NanoPlot --listcolors
NanoPlot --listcolormaps
NanoPlot --version

# Merged tests - combining multiple options per execution:

echo "testing bam with pdf format, outlier dropping, title, and quality filtering:"
NanoPlot --bam nanotest/alignment.bam --verbose -o tests --format pdf --drop_outliers --title "Test Run" --minqual 8

echo "testing fasta with tsv_stats, length filtering, and huge mode:"
NanoPlot --fasta nanotest/reads.fa.gz --verbose --maxlength 35000 -o tests --tsv_stats --huge

echo "testing bam without supplementary, with colormap and no_static:"
NanoPlot --bam nanotest/alignment.bam --verbose --no_supplementary -o tests --colormap Viridis --no_static

echo "testing summary with loglength, N50, multiple plots, and info_in_report:"
NanoPlot --summary nanotest/sequencing_summary.txt --loglength --verbose -o tests --N50 --plots kde dot --info_in_report

echo "testing fastq_rich with downsample, store, and custom prefix:"
NanoPlot --fastq_rich nanotest/reads.fastq.gz --verbose --downsample 500 -o tests --store --prefix test_

echo "testing fastq_minimal with plots, color, threads, and only-report:"
NanoPlot --fastq_minimal nanotest/reads.fastq.gz --verbose --plots dot -o tests --color red --threads 2  --only-report

# Check if CRAM file exists before testing
if [ -f "nanotest/alignment.cram" ]; then
    echo "testing cram with multiple formats, length filtering:"
    NanoPlot --cram nanotest/alignment.cram --verbose -o tests --format png jpeg --minlength 1000 --maxlength 40000
else
    echo "CRAM file not found, skipping CRAM test (run make_cram.sh to create it)"
fi

# Optional legacy test (commented out due to dependency issues)
# echo "testing legacy with summary:"
# pip install seaborn==0.10.1 "numpy<1.24"
# NanoPlot --summary nanotest/sequencing_summary.txt --loglength --verbose -o tests --legacy hex --raw --prefix legacy_ --plots dot --format pdf