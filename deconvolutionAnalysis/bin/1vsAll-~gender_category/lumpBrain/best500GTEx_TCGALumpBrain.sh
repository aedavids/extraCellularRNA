#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 11/15/23
# 
# pipeline.sh is generic/reuseable.
# to reduce user and environment specifc coupling it takes a lot of arguments
# This script make calling the pipeline easier and makes the results reproducable
# we capture all the parameter values
#
# ref: extraCellularRNA/deconvolutionAnalysis/bin/test.slurm.sh
#
#

if [ $# -ne 2 -a $# -ne 3 ];
    then
        printf "ERROR  \n"
        printf "local install: cp ~/extraCellularRNA/deconvolutionAnalysis/bin/{1vsAll-~gender_category/${0},pipeline.sh} .\n"
        printf "usage:"
        printf "to run on current server: $0 ciberSortSecurityToken ciberSortUser\n"
        printf "to run on slurm         : $0 ciberSortSecurityToken ciberSortUser SLURM\n"
        printf "usage: tail -f ${0}.log \n"
        printf "follow 'Token and instruction access' @ https://cibersortx.stanford.edu/download.php \n"
        exit 1 # error
    fi

#
# Copy and past: Common paramters to change
# - topN
#


ciberSortSecurityToken=$1
ciberSortUser=$2
runSlurm=$3

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

dataSet='GTEx_TCGA'

# rootDir="/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python"
rootDir="/private/groups/kimlab/GTEx_TCGA"

# colData="${rootDir}/pipeline/dataFactory/test/data/testIntegration/colData.csv"
# colData="${rootDir}/groupbyGeneTrainingSets/GTEx_TCGA_TrainColData.csv"
colData="${rootDir}/groupbyGeneTrainingSets/${dataSet}_TrainLumpBrain.colData.csv"

# countData="${rootDir}/pipeline/dataFactory/test/data/testIntegration/geneCounts.csv"
countData="${rootDir}/groupbyGeneTrainingSets/GTEx_TCGA_TrainGroupby.csv"

# deseqResultsDir="${rootDir}/pipeline/dataFactory/test/data/testSignatureGenes/1vsAll" 
# ~gender + category 
# deseqResultsDir="${rootDir}/1vsAll"
deseqResultsDir="${rootDir}/1vsAllLumpBrain"

# findModule="analysis.createBest25CreateSignatureGeneConfig"
findModule="analysis.createBestCreateSignatureGeneConfig"

# estimatedScalingFactors="${rootDir}/pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/estimatedSizeFactors.csv"
estimatedScalingFactors="${rootDir}/1vsAll/estimatedSizeFactors.csv"

outDir=`pwd`"/${0}.out"

# clean up any leftovers from old runs
'rm' -f serial*.log;

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
#
wdlRoot="${gitRoot}/terra/cibersortx/wdl"


rm -rf wdlTest
mkdir wdlTest
'cp' "${wdlRoot}/CIBERSORTxFractionsWorkflow.wdl" .
# -j just store file name ignore directory path
zip -j import.zip "${wdlRoot}/cibersortxFractionsTask.wdl" "${wdlRoot}/../wdlTest/partitionDataTask.wdl" "${wdlRoot}/../wdlTest/mergeTask.wdl"

wdlInputJSON="${outDir}/CIBERSORTxFractionsWorkflow.wdl.input.json"

# arguments to pass to BestSignatureGeneConfig.__init__()
topN=500
title=` printf "best%s" $topN` # do not uses spaces, title will be part of file paths
vargs=" --design tilda_gender_category \
        --padjThreshold 0.001 \
        --lfcThreshold 2.0 \
        --dataSetName  $dataSet\
        --number  $topN \
        --title ${title} \
        --localCacheRoot ${outDir} "

printf "\n\n\n SignatureGeneConfig vargs : $vargs \n !!!!!! \n"

logFile="${0}.log"
rm -f $logFile
if [ -z "${runSlurm}" ]; then
    # runSlurm is empty, start pipeline as background process
    setsid sh -c "set -x; pipeline.sh ${colData} \
                            ${countData} \
                            ${deseqResultsDir} \
                            ${findModule} \
                            ${estimatedScalingFactors} \
                            ${outDir} \
                            import.zip \
                            ${wdlInputJSON} \
                            ${WDL_TOOLS} \
                            ${gitRoot} \
                            ${ciberSortSecurityToken} \
                            ${ciberSortUser} \
                            ${vargs}" > $logFile 2>&1 & 

    sleep 10
    pstree $USER
    ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep $USER
else
    if [[ $runSlurm != "SLURM" ]]; then
        printf "ERROR the 3rd argument $runSlurm != SLURM"
    else
        sbatch pipeline.sh ${colData} \
                            ${countData} \
                            ${deseqResultsDir} \
                            ${findModule} \
                            ${estimatedScalingFactors} \
                            ${outDir} \
                            import.zip \
                            ${wdlInputJSON} \
                            ${WDL_TOOLS} \
                            ${gitRoot} \
                            ${ciberSortSecurityToken} \
                            ${ciberSortUser} \
                            ${vargs} > $logFile 2>&1
    fi
fi