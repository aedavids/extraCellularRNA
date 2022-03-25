#!/bin/bash

scriptName=`basename $0`
if [ $# -ne 2 ] ; then
    echo "ERROR: usage $scriptName destinationBucketId  quantFile.csv"
    echo "example: $scriptName  gs://anvil_gtex_v8_hg38_edu_ucsc_kim_lab_spark terraDataModels/test-aedavids-proj/AnVIL_GTEx_V8_hg38_edu_ucsc_kim_lab/createDESeqData.output/GTEx-trainQuantFiles.csv"
    echo "quantFile.csv has two columns 'sampleName,source'"
    exit 1
fi

#set -x # turn debug on
# set +x # turn debug off

bucketURL="$1"
colDataFile="$2"

echo "AEDWIP remove head"

# use sed to remove the first line. it is a header
data=`sed  1d  < $colDataFile`

for i in $data;
do
    sampleId=`echo $i |cut -d , -f 1`
    url=`echo $i |cut -d , -f 2`

    set -x
    gsutil -m cp $url gs://${bucketURL}/quant/${sampleId}.quant.sf.gz
    set +x 
done

