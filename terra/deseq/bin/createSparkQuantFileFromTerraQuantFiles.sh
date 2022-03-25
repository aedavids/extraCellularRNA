#!/bin/bash

#
# the urls in the
# terraDataModels/test-aedavids-proj/AnVIL_GTEx_V8_hg38_edu_ucsc_kim_lab/createDESeqData.output/
# quantFiles.csv point to our terra workspace bucket. Spark must run in a native GCP project
# and can not access these file. We used copyQuantFilesFromTerra.sh to copy the to the native GCP project
#

scriptName=`basename $0`
if [ $# -ne 2 ] ; then
    echo "ERROR: usage $scriptName  bucketId quantFile.csv"
    echo "example: $scriptName anvil-gtex-v8-hg38-edu-ucsc-kim-lab-spark terraDataModels/test-aedavids-proj/AnVIL_GTEx_V8_hg38_edu_ucsc_kim_lab/createDESeqData.output/GTEx-trainQuantFiles.csv"
    echo "quantFile.csv has two columns 'sampleName,source'"
    exit 1
fi

#set -x # turn debug on
# set +x # turn debug off

bucket="$1"
colDataFile="$2"

header=`head -n1 $colDataFile`
echo $header

# use sed to remove the first line. it is a header
data=`sed  1d  < $colDataFile `


for i in $data;
do
    sampleId=`echo $i |cut -d , -f 1`
    url=`echo $i |cut -d , -f 2`

    #set -x
    #gsutil -m cp $url gs://${bucketURL}/quant/${sampleId}.quant.sf.gz
    echo "${sampleId},gs://${bucket}/quant/${sampleId}.quant.sf.gz"
    #set +x 
done

