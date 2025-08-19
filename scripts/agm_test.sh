#! /bin/bash

#SBATCH --time=04-00:00:00
#SBATCH --partition=defq
#SBATCH --mail-user=myemail@email.org
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks-per-node=64
#SBATCH --mem=128GB
#SBATCH --nodes=1
#SBATCH --job-name=nplot
#SBATCH --comment=nplot

source /path/to/nanoplot_env/bin/activate

# test fresh nanoplot with update
python /path/to/NanoPlot/nanoplot/NanoPlot.py --fastq /path/to/test_file.fastq.gz --verbose --minqual 4 --color red -o scripts/agm_tests
