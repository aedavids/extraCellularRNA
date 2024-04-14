#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 4/1/24
# 
# runs  extraCellularRNA/intraExtraRNA_POC/python/src/models/randomForestHyperparmeter 
# on elife Esophagus cancer samples
#

# 
# to use slurm see extraCellularRNA/deconvolutionAnalysis/bin/test.slurm.sh
#
# ref:
#   https://giwiki.gi.ucsc.edu/index.php/Overview_of_using_Slurm
#   extraCellularRNA/terra/cibersortx/wdl/README.md 
#   extraCellularRNA/terra/cibersortx/wdl/CIBERSORTxFractionsWorkflow.slurm.sh
#
# slurm argument are ignored.
#
#
# Partition - This is the queue it goes in:
#SBATCH --partition=medium
#
# Where to send email (optional)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aedavids@ucsc.edu
#

#
# configure required resoures
#
# Number of nodes you need per job:
#SBATCH --nodes=1
#
# Memory needed for the jobs.  Try very hard to make this accurate.  DEFAULT = 4gb
# Default units are megabytes. Different units can be specified using the suffix [K|M|G|T
# S-aedwip-BATCH --mem=4gb
#S-aedwip-OOM-BATCH --mem=32G # total memory per node
#SBATCH --mem=18M # total memory per node
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=1
#
# Processors per task:
# At least eight times the number of GPUs needed for nVidia RTX A5500
# ref https://giwiki.gi.ucsc.edu/index.php?title=Using_Docker_under_Slurm
# slurm cannot limit the resources that docker uses. 
# set to number of samples / CIBERSORTxFractionsWorkflow.wdl.input.json CIBERSORTxFractionsWorkflow.numSamplesInPartition
#SBATCH --cpus-per-task=32
#
# Number of GPUs, this can be in the format of "--gres=gpu:[1-8]", or "--gres=gpu:A5500:[1-8]" with the type included (optional)
#SBATCH --gres=gpu:1
#
# Standard output and error log
#SBATCH --output=serial_test_%j.log 
#
# Wall clock limit in hrs:min:sec:
# best500GTEx* timmed out
#S-AEDWIP-BATCH --time=10:02:30
#SBATCH --time=24:00:00
#

#
## Command(s) to run (example):
#


if [ $# -ne 1 -a $# -ne 2 ];
    then
        printf "ERROR  \n"
        printf "local install: cp ~/extraCellularRNA/intraExtraRNA_POC/bin/$0 .\n"
        printf "usage to run on current server:  $0 outDir \n"
        printf "usage to run on SLURM         :  $0 outDir SLURM\n"
        exit 1 # error
    fi

progName=`basename $0`
outDir=$1/${progName}.out
mkdir $outDir

runSlurm=$2

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options


#
# set up the python env
#
printf "\n\n\nset up python environment \n"
# start conda env
condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`
source ${condaBase}/etc/profile.d/conda.sh
conda activate extraCellularRNA

pythonSrcRoot="/private/home/aedavids/extraCellularRNA"

if [ -z ${PYTHONPATH+x} ];
    then
        #PYTHONPATH is unset or set to the empty string:
        export PYTHONPATH="${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
        export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/intraExtraRNA_POC/python/src"; 
    else 
        export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
        export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/intraExtraRNA_POC/python/src"; 
    fi

printf "PYTHONPATH : $PYTHONPATH \n"

logFile="${outDir}/${progName}.log"

if [ -z "${runSlurm}" ]; then
    # runSlurm is empty, start pipeline as background process
    setsid sh -c "set -x; python -m models.randomForestHyperparmeterSearch \
                        --outDir ${outDir} \
                        --features Esophagus_Muscularis \
                        --elife 'Healthy donor' 'Esophagus Cancer'" > $logFile 2>&1 & 
else
    if [[ $runSlurm != "SLURM" ]]; then
        printf "ERROR the 2nd argument $runSlurm != SLURM"
    else
        todo this does not work. we need activate our conda env
        sbatch --job-name ${progName} sh -c "set -x; python -m models.randomForestHyperparmeterSearch \
                    --outDir ${outDir} \
                    --features Esophagus_Muscularis \
                    --elife 'Healthy donor' 'Esophagus Cancer'" > $logFile 2>&1 &
    fi
fi

sleep 10
pstree $USER
 ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep $USER
