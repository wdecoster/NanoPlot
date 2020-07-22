set -ev

git clone https://github.com/wdecoster/nanotest.git

NanoPlot -h
NanoPlot --listcolors
echo ""
echo ""
echo ""
echo "testing bam:"
NanoPlot --bam nanotest/alignment.bam --verbose
echo ""
echo ""
echo ""
echo "testing bam without supplementary alignments:"
NanoPlot --bam nanotest/alignment.bam --verbose --no_supplementary
echo ""
echo ""
echo ""
echo "testing summary:"
NanoPlot --summary nanotest/sequencing_summary.txt --loglength --verbose
echo ""
echo ""
echo ""
echo "testing fastq rich:"
NanoPlot --fastq_rich nanotest/reads.fastq.gz --verbose --downsample 800
echo ""
echo ""
echo ""
echo "testing fastq minimal:"
NanoPlot --fastq_minimal nanotest/reads.fastq.gz --store --verbose --plot dot
echo ""
echo ""
echo ""
echo "testing fastq plain:"
NanoPlot --fastq nanotest/reads.fastq.gz --verbose --minqual 4 --color red
echo ""
echo ""
echo ""
echo "testing fasta:"
NanoPlot --fasta nanotest/reads.fa.gz --verbose --maxlength 35000
echo ""
echo ""
echo ""
echo "testing feather:"
NanoPlot --feather nanotest/summary1.feather --verbose --outdir plots
