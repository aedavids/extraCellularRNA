#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 1/14/2025
# 

#
# output env info to make debugging easier
#
set -x
scriptName=`basename $0`
pwd; 
hostname; 
date

set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

img="aedavids/edu_ucsc_kim_lab-1vsall_1.0"

outputDir="`pwd`/${scriptName}.output"
mkdir -p "${outputDir}"
printf "\n\n\n************ AEDWIP do not run docker we do not have a valid token\n"

testDataRoot="/private/home/aedavids/extraCellularRNA/terra/deseq/R/data/1vsAllTest"

#
# you can debug docker token problems by removing
# the --detach flag, --rm and adding --interactive --tty
# this will cause container error message to be written to the terminal
#
USER_ID=`id -u`

#         --design  '\"~ sex + tissue_id\"'   


#
# the following debug from cli works
# start docker
# run deseq on docker
# testDataRoot="/private/home/aedavids/extraCellularRNA/terra/deseq/R/data/1vsAllTest"

# outputDir="/scratch/aedavids/aedwip.output"
# mkdir -p $outputDir
# /home/rstudio/DESeqScript.R         --countMatrix /data/unitTestGroupByGenesCountMatrix.csv         --colData /data/unitTestColData.csv         --design  '~ sex + tissue_id'         --referenceLevel Lung         --outFile /outDir/aedwipResults.tsv         --estimateSizeFactorsOutfile /outDir/aedwipScalingFactors.tsv         --isCSV
# docker run  --interactive --tty -e USERID=30108 -v ${testDataRoot}:/data -v ${outputDir}:/outDir aedavids/edu_ucsc_kim_lab-1vsall_1.0 /bin/bash
cmd="docker run \
    --interactive --tty \
    -e USERID=${USER_ID} \
    -v ${testDataRoot}:/data \
    -v ${outputDir}:/outDir \
    ${img} \
    /home/rstudio/DESeqScript.R \
        --countMatrix /data/unitTestGroupByGenesCountMatrix.csv \
        --colData /data/unitTestColData.csv \
        --referenceLevel Lung \
        --outFile /outDir/aedwipResults.tsv \
        --estimateSizeFactorsOutfile /outDir/aedwipScalingFactors.tsv \
        --isCSV \ 
        --design aedwipDesign
    "

#    --detach \
#     --rm \

    # -v ${cibersortInputDir}:/src/data \
    # -v ${cibersortOutputDir}:/src/outdir \
    # ${img} \
    # --username ${cibersortUser} \
    # --token ${cibersortToken} \
    # --mixture ${mixtureMatrix}\
    # --sigmatrix ${signatureMatrix}\
    # --perm 100 \
    # --label $jobId \
    # --QN FALSE \
    # --verbose TRUE

printf "\n\n\n************ run docker\n"
echo "${cmd}" > "${scriptName}.docker.parameters.txt"

# scriptOut="$cibersortInputDir/${scriptName}.meta.out"
# echo "run on ${timeStamp}" > ${scriptOut}
# echo "input src: ${bestSrc}/ciberSort/*"  >> ${scriptOut}
# echo $cmd >> ${scriptOut}
# echo ""   >> ${scriptOut}

echo ""
$cmd 
exitStatus=$?
if [ $exitStatus -ne 0 ]; then
    printf "docker run failed with exit status : ${exitStatus}\n"
    exit $exitStatus
fi
