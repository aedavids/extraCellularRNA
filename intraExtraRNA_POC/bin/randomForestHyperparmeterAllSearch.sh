#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 2/5/24
# 
# runs  extraCellularRNA/intraExtraRNA_POC/python/src/models/randomForestHyperparmeter 
# on elife lung cancer samples
#


if [ $# -ne 1 ];
    then
        printf "ERROR  \n"
        printf "local install: cp ~/extraCellularRNA/intraExtraRNA_POC/bin/$0 .\n"
        printf "usage: $0 outDir \n"
        exit 1 # error
    fi

progName=`basename $0`
outDir=$1/${progName}.out
mkdir $outDir

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
setsid sh -c "set -x; python -m models.randomForestHyperparmeterSearch \
                        --outDir ${outDir} \
                        --features Lung LUAD LUSC \
                        --elife 'Healthy donor' 'Lung Cancer' \
                        " > $logFile 2>&1 & 

sleep 10
pstree $USER
 ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep $USER
