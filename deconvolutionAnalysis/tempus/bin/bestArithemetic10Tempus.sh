#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 12/28/2024
# 
# runs deconvolutionPipeline.sh
# selects DESeq2 genes based on abs(LFC) > 2.0 , padj=0.001. It sorts 
# the results by base mean in descending order and selects the top 10.
#
# deconvolutionPipeline.sh is generic/reuseable.
# to reduce user and environment specifc coupling it takes a lot of arguments
# This script make calling the pipeline easier and makes the results reproducable
# we capture all the parameter values
#
#

scriptName=`basename $0`
if [ $# -ne 3 ];
    then
        printf "ERROR  \n"
        printf "local install: cp ~/extraCellularRNA/deconvolutionAnalysis/python/tempus/bin/{${scriptName},deconvolutionArithemeticPipeline.sh} .\n"
        printf "usage: $0 ciberSortUser ciberSortSecurityToken tumorId \n"
        printf "usage: tail -f ${0}.log \n"
        printf "follow 'Token and instruction access' @ https://cibersortx.stanford.edu/download.php \n"
        exit 1 # error
    fi

#
# Copy and past: Common paramters to change
# - topN
#

set -x
topN=10
ciberSortSecurityToken=$1
ciberSortUser=$2
tumorId=$3
#
# cibersort input and output directory mount points
# You must use full paths. For unknown reasons cibersort raises an error if 
# the output directory in my home directory.
# the output and input directories can be the same
cibersortInputDir="${PWD}/${scriptName}.out/cibersortInputDir"
# cibersortOutputDir=/scratch/aedavids/cibersortOut 
cibersortOutputDir="${PWD}/${scriptName}.out/cibersortOut" 
mkdir -p "${cibersortOutputDir}" "${cibersortOutputDir}"


normalizedCountFilePath="/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv"
metaDataFilePath="/private/groups/kimlab/data/tempus/illumina/20241107/raw/metaDataWithHeader.csv"

#
######################################################################
# run the pipeline


# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

outDir=${scriptName}.out
mkdir -p $outDir

logFile="${0}.log"
rm -f $logFile
setsid sh -c "set -x; deconvolutionArithemeticPipeline.sh \
    ${ciberSortSecurityToken} \
    ${ciberSortUser} \
    ${topN} \
    ${normalizedCountFilePath} \
    ${metaDataFilePath}  \
    ${tumorId} \
    ${cibersortInputDir} \
    ${cibersortOutputDir}" \
        > $logFile 2>&1 & 

#
# display pipeline process information
# it may take a couple of seconds for the docker container to start
sleep 10
pstree $USER
ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep $USER

