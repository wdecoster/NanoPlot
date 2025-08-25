#! /bin/bash

#SBATCH --time=04-00:00:00
#SBATCH --partition=defq
#SBATCH --mail-user=email@email.org
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks-per-node=64
#SBATCH --mem=128GB
#SBATCH --nodes=1
#SBATCH --job-name=nplot
#SBATCH --comment=nplot

source /home/tmhagm8/scratch/nanoplot_env/bin/activate

# Go to the repo root
cd /home/tmhagm8/scratch/NanoPlot

# Make sure to use right Python imports
export PYTHONPATH="$PWD:$PYTHONPATH"

# Double check imports
python - <<'PY'
import nanoplotter.plot as p
import nanoplot.utils as u
print("USING nanoplotter.plot:", p.__file__)
print("USING nanoplot.utils :", u.__file__)
PY

# check it 
python -m nanoplot.NanoPlot \
  --fastq /home/tmhagm8/scratch/SOMAteM_bckp/SOMAteM/examples/data/B011_2.fastq.gz \
  -t 14 --verbose --minqual 4 --dpi 600  \
  -o /home/tmhagm8/scratch/NanoPlot/scripts/agm_test -f png
