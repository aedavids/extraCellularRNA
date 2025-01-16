#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 1/15/2025
# 
# assumes each tummor sample is unique
# 1. select samples for tumor id and create count matrix we can use with DESeq2
# 2. run DESeq2
#   we have hack the count and meta data to get DESeq2 to run
#
#   a) work around for error 
#       'The design matrix has the same number of samples and coefficients to fit ...'
#   we only have 2 samples and DESeq2 . 
#   we add a fake control and a fake UD sample
#
#   b) work around for error 
#   'Error in estimateDispersionsFit(object, fitType = fitType, quiet = quiet) : 
#   all gene-wise dispersion estimates are within 2 orders of magnitude ...'
#   we add error handling to the DESeq2 script. Will use gene-wise dispersion 
#   estimates if the standard dispersion fitting fails


#
# output env info to make debugging easier
#
scriptName=`basename $0`
pwd; 
hostname; 
date

set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

#
# parse arguments
#

# print all the cli argumnents
printf "################ ${scriptName} : command line arguments \n"
for var in "$@"
do
    printf "argument :$var \n"
done

# parse the arguments
numberOfArguments=1

if [ $# -ne $numberOfArguments ];
    then
        printf "ERROR ${scriptName} missing command line arguments. expected $numberOfArguments recevied $# \n"
        exit 1 # error
    fi

tokenId=$1 # example T1 

# outDir="${scriptName}.out"
# mkdir -p "${outDir}"
# pushd "${outDir}"

#
################## step 1: create the count matrix
# 
printf "\n\n\n################ step 1: create the count matrix \n"
createTumorCountMatrix.sh
countDataDir=createTumorCountMatrix.sh.output

#
################## step 2: run DESeq2
#
printf "\n\n\n################ step 2: run DESeq2 \n"
img="aedavids/edu_ucsc_kim_lab-1vsall_1.1"

# docker can not write to /private/home/aedavids
#deseqOutDir="edu_ucsc_kim_lab-1vsall_1.1.out"
deseqOutDir="/scratch/aedavids/edu_ucsc_kim_lab-1vsall_1.1.out"

mkdir -p "${deseqOutDir}"

USER_ID=`id -u`
design="~ category"

DESEq2InputDir="${PWD}/DESeq2Input"
mkdir -p "${DESEq2InputDir}"

metaDataFilePath="/private/groups/kimlab/data/tempus/illumina/20241107/raw/metaDataWithHeader.csv"
head -n 1 "${metaDataFilePath}" > "${DESEq2InputDir}/metaDataWithHeader.csv"
grepPatern="${tokenId}_UD\|${tokenId}_Control"
printf "grepPatern: ${grepPatern}\n"
grep "${grepPatern}" "${metaDataFilePath}" >> "${DESEq2InputDir}/metaDataWithHeader.csv"

# 'cp' "${metaDataFilePath}" "${DESEq2InputDir}"
chmod a+w "${DESEq2InputDir}/metaDataWithHeader.csv"

DESeq2CountMatrix="${tokenId}RawDESeq2FmtControlUDCounts.csv"
'cp' "${countDataDir}/${DESeq2CountMatrix}" "${DESEq2InputDir}"

# DESeq2 throws error.
# this is because we only 2 samples one is the control the other is UD
# and our design has 2 parameter
#
# Error in checkForExperimentalReplicates(object, modelMatrix) : 
#
#   The design matrix has the same number of samples and coefficients to fit,
#   so estimation of dispersion is not possible. Treating samples
#   as replicates was deprecated in v1.20 and no longer supported since v1.22.
#
# Calls: estimateDispersions ... estimateDispersions -> .local -> checkForExperimentalReplicates
#
# the following hack adds a fake control and a fake UD sample
#

pushd "${DESEq2InputDir}"

printf "\nadd fake samples to meta data\n" > t
tail -n 2 metaDataWithHeader.csv | sed -e s/SLD/fake/g > t
cat metaDataWithHeader.csv t > tt
'mv' tt metaDataWithHeader.csv

printf "\nadd fake samples to count data\n" > t
cut -d , -f 2,3 "${DESeq2CountMatrix}"  | sed -e 's/SLD/fake/g' > t; 
paste -d ,  "${DESeq2CountMatrix}" t > tt; 
'mv' tt "${DESeq2CountMatrix}"

popd

#
# you can debug docker token problems by removing
# the --detach flag, --rm and adding --interactive --tty
# this will cause container error message to be written to the terminal
#

    # --detach \
    # --rm \

cmd="docker run \
    --interactive --tty \
    -e USERID=${USER_ID} \
    -v ${DESEq2InputDir}:/data \
    -v ${deseqOutDir}:/outDir \
    ${img} \
    /home/rstudio/DESeqScript.R \
        --countMatrix /data/${DESeq2CountMatrix} \
        --colData /data/metaDataWithHeader.csv \
        --referenceLevel Control \
        --outFile /outDir/${tokenId}Results.tsv \
        --estimateSizeFactorsOutfile /outDir/${tokenId}ScalingFactors.tsv \
        --isCSV \
        --design \"${design}\"
    "

#
# in all are other docker scripts we format a string then use
# $cmd to run it
# because design is a variable argument list this does not work
# we can not pass the design string correctly
# we get the error '/home/rstudio/DESeqScript.R: error: unrecognized arguments: sex + tissue_id"
#
# the eval command is needed to expand the quotes in the cmd string
# $cmd
eval $cmd 
exitStatus=$?
if [ $exitStatus -ne 0 ]; then
    printf "docker run failed with exit status : ${exitStatus}\n"
    exit $exitStatus
fi

printf "\n\n################ ${scriptName} END \n"
printf "PWD: $PWD \n"
printf "checkoutput directories: ${countDataDir} ${DESEq2InputDir} ${deseqOutDir} \n"