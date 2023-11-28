#!/bin/bash
# https://giwiki.gi.ucsc.edu/index.php/Overview_of_using_Slurm
#
# Job name:
#SBATCH --job-name=aedavids-wdlTest
#
# Partition - This is the queue it goes in:
#SBATCH --partition=main
#
# Where to send email (optional)
#SBATCH --mail-user=aedavids@ucsc.edu
#
# Number of nodes you need per job:
#SBATCH --nodes=1
#
# Memory needed for the jobs.  Try very hard to make this accurate.  DEFAULT = 4gb
#SBATCH --mem=32gb
#
# Number of tasks (one for each CPU desired for use case) (example):
# https://stackoverflow.com/a/53759961
# The --ntasks parameter is useful if you have commands that you want to run in parallel
# within the same batch script. This may be two separate commands separated by an & or
# two commands used in a bash pipe (|).
#SBATCH --ntasks=1
#
# Processors per task:
# At least eight times the number of GPUs needed for nVidia RTX A5500
#SBATCH --cpus-per-task=32
#
# Number of GPUs, this can be in the format of "--gres=gpu:[1-8]", or "--gres=gpu:A5500:[1-8]" with the type included (optional)
##SBATCH --gres=gpu:1
#
# Standard output and error log
#SBATCH --output=serial_test_%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=06:00:30
#
## Command(s) to run (example):
pwd; hostname; date

# start conda env
condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`
source ${condaBase}/etc/profile.d/conda.sh
set -x
conda activate extraCellularRNA

which java

java --version

docker pull aedavids/wdltest
docker pull aedavids/cibersortx_fractions

runCIBERSORTxFractionsWorkflow.sh
echo "aedavids-CIBERSORTxFractionsWorkflow.slurm.sh done!"
date
