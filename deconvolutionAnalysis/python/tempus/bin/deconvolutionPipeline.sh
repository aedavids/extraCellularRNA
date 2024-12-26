#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 12/25/2024
# 

# TODO AEDWIP parse command line arguments
scriptName=`basename $0`

lfcThreshold=2.0 
topN=10
padjThreshold=0.001
deseq2ResultsFilePath="/private/groups/kimlab/data/tempus/illumina/20241107/create/results/data/Control_vs_UD_deseq_results.csv"
outDir="./tmp"


set -x turn debug trace on. output goes to stderr, normal output goes to stdout

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
# set -x
conda activate extraCellularRNA

pythonSrcRoot="/private/home/aedavids/extraCellularRNA"

if [ -z ${PYTHONPATH+x} ];
    then
        #PYTHONPATH is unset or set to the empty string:
        export PYTHONPATH="${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
    else 
        export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
    fi

printf "PYTHONPATH : $PYTHONPATH \n"


#
# find biomarkers
# 
python -m tempus.findBiomarkers  \
    --lfcThreshold ${lfcThreshold} \
    --number ${topN} \
    --padjThreshold ${padjThreshold} \
    --deseq2ResultsFilePath ${deseq2ResultsFilePath} \
    --outDir ${outDir} 

#
# get the top N biomarkers
# use cut to get the gene_id column
# grep -v remove the column header
# xargs echo -n remove the newline
genesOfInterest=`cat ${outDir}/biomarkerDESeq2Results.csv | cut -d , -f 7 | grep -v gene_id | xargs echo -n`
printf "genesOfInterest : ${genesOfInterest}\n"
