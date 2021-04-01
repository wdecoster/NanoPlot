set -ev

if [ -d "nanotest" ]; then
    echo "nanotest already cloned"
else
    git clone https://github.com/wdecoster/nanotest.git
fi

NanoPlot -h
# NanoPlot --listcolors
echo ""
echo ""
echo ""
echo "testing bam:"
NanoPlot --bam nanotest/alignment.bam --verbose -o tests
echo ""
echo ""
echo ""
echo "testing bam without supplementary alignments:"
NanoPlot --bam nanotest/alignment.bam --verbose --no_supplementary -o tests
echo ""
echo ""
echo ""
echo "testing summary:"
NanoPlot --summary nanotest/sequencing_summary.txt --loglength --verbose -o tests
echo ""
echo ""
echo ""
echo "testing fastq rich:"
NanoPlot --fastq_rich nanotest/reads.fastq.gz --verbose --downsample 800 -o tests
echo ""
echo ""
echo ""
echo "testing fastq minimal:"
NanoPlot --fastq_minimal nanotest/reads.fastq.gz --store --verbose --plot dot -o tests
echo ""
echo ""
echo ""
echo "testing fastq plain:"
NanoPlot --fastq nanotest/reads.fastq.gz --verbose --minqual 4 --color red -o tests
echo ""
echo ""
echo ""
echo "testing fasta:"
NanoPlot --fasta nanotest/reads.fa.gz --verbose --maxlength 35000 -o tests
echo ""
echo ""
echo ""
# echo "testing feather:"
# NanoPlot --feather nanotest/summary1.feather --verbose --outdir plots
