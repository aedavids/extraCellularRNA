#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 11/15/23
# 
# pipeline.slurm.sh is generic/reuseable.
# to reduce user and environment specifc coupling it takes a lot of arguments
# This script make calling the pipeline easier and makes the results reproducable
# we capture all the parameter values
#

#
# copy everything to phoenix
#  cp ~/extraCellularRNA/deconvolutionAnalysis/bin/{test.slurm.sh,pipeline.slurm.sh} .; ./test.slurm.sh aedwipSecurityToken ${USER}@ucsc.edu
#

if [ $# -ne 2 ];
    then
        printf "ERROR  \n"
        printf "usage: $0 ciberSortSecurityToken ciberSortUser\n"
        printf "follow 'Token and instruction access' @ https://cibersortx.stanford.edu/download.php"
        exit 1 # error
    fi

ciberSortSecurityToken=$1
ciberSortUser=$2

rootDir="/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python"
colData="${rootDir}/pipeline/dataFactory/test/data/testIntegration/colData.csv"
countData="${rootDir}/pipeline/dataFactory/test/data/testIntegration/geneCounts.csv"
deseqResultsDir="${rootDir}/pipeline/dataFactory/test/data/testSignatureGenes/1vsAll" 
findModule="pipeline.dataFactory.test.exampleCreateSignatureGeneConfig"
estimatedScalingFactors="${rootDir}/pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/estimatedSizeFactors.csv"
outDir=`pwd`"/pipeline.slurm.sh.out"

# clean up any leftovers from old runs
'rm' -i serial*.log;

printf "deleted ${outDir}?\n"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) 'rm' -rf "${outDir}"; break;;
        No )  break;;
    esac
done

printf "deleted cromwell output?\n"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) 'rm' -rf cromwell-* ; break;;
        No )  break;;
    esac
done

'rm' -f import.zip
'rm' -f CIBERSORTxFractionsWorkflow.wdl
'rm' -f cromwellOptions.json 

mkdir -p "${outDir}"

printf "\n\n\n"

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options


#
# prepare cromwell/wdl workflow input
# The pipeline will call cromwell and pass the cibersort input files
#
gitRoot="/private/home/aedavids/extraCellularRNA"
WDL_TOOLS="${gitRoot}/java/bin"


'cp' "${gitRoot}/deconvolutionAnalysis/bin/cromwellOptions.json.template" .

#
# copy required wdl files 
# AEDWIP TODO check if we need zip file or not?
# wdl import are ../wdlTest/my.wdl
# should we change to import my.wdl ?
#
wdlRoot="${gitRoot}/terra/cibersortx/wdl"


rm -rf wdlTest
mkdir wdlTest
'cp' "${wdlRoot}/CIBERSORTxFractionsWorkflow.wdl" .
# -j just store file name ignore directory path
zip -j import.zip "${wdlRoot}/cibersortxFractionsTask.wdl" "${wdlRoot}/../wdlTest/partitionDataTask.wdl" "${wdlRoot}/../wdlTest/mergeTask.wdl"

wdlInputJSON="${outDir}/CIBERSORTxFractionsWorkflow.wdl.input.json"

sbatch pipeline.slurm.sh "${colData}" \
                        "${countData}" \
                        "${deseqResultsDir}" \
                        "${findModule}" \
                        "${estimatedScalingFactors}" \
                        "${outDir}" \
                        import.zip \
                        "${wdlInputJSON}" \
                        "${WDL_TOOLS}" \
                        "${gitRoot}" \
                        "${ciberSortSecurityToken}" \
                        "${ciberSortUser}"

squeue
sleep 2

# if tail fails it is probably because the job is in the queue 
# and has not started running yet
tail -f serial*.log
