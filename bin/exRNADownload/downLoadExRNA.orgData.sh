#!/bin/bash

# download data from exRNA.org
# aedavids@ucsc.edu

export scriptName=`basename $0`

if [ $# -ne 6 ]; then
    echo "error: usage $scriptName rootDir downLoadFile url sampleId condition assession"
    echo "example $scriptName ./tmp Sample_4S2.fastq.zip ftp://ftp.genboree.org/exRNA-atlas/grp/Extracellular%20RNA%20Atlas/db/exRNA%20Repository%20-%20hg19/file/exRNA-atlas/exceRptPipeline_v4.6.2/TPATE1-human-plasma-healthyVsCancer-2016-10-17/sample_Sample_4S2_fastq/rawInput/Sample_4S2.fastq.zip Sample_4S2 Colon_Carcinoma EXR-TPATE1NEBLIB4S2-BS"
    exit 1
fi

rootDir=$1
downloadFileName=$2
url=$3
sampleId=$4
condition="$5"
assession="$6"

#echo ""
#echo "rootDir           $rootDir"
#echo "downloadFileName  $downloadFileName"
#echo "url               $url"
#echo "sampleId          $sampleId"
#echo "condition         $condition"
#echo "assession         $assession"

downLoadDir="${rootDir}/${condition}/${sampleId}"
mkdir -p "${downLoadDir}"

dataFile="${downLoadDir}/${downloadFileName}"
if test -f "$dataFile"; then
    echo "WARNING $dataFile exists. download was skipped."
    exit 0
fi

echo ""
set -x # turn debug trace on
# --silent do not show progress. prevent curl from generating huge log files
curl $url --silent --output "${dataFile}"
existStatus=$?
set +x # turn debug trace off
echo ""

exit $exitStatus
