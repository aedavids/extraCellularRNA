#!/bin/bash

# aedavids@ucsc.edu
## download tools and install
# - [SRA toolkit home page](https://ncbi.github.io/sra-tools/)
# - [download tool for centos](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz)
# - [how to use sra tools](https://edwards.sdsu.edu/research/fastq-dump/)


export scriptName=`basename $0`
if [ $# -ne 2 ]; then
    echo "error: usage $scriptName phoneNumber fileOfSamples"
    echo "example $scriptName 7508622639@txt.att.net SRR_Acc_List.txt.part.ad"
    echo "SRR_Acc_List.txt.part.ad is a list of samples"
    echo "fileOfSample creation example: 'split -l 100 SRR_Acc_List.txt SRR_Acc_List.txt.part.' "
    exit 1
fi


# use export else variable will not be defined by inline script embedded in nohup
export phoneNumber=$1
export fileOfSamples=$2
export listOfSamples=`cat ${fileOfSamples}`
export dateStamp=`dateStamp.sh`

set -x
which fastQDumpRunner.sh
nohup sh -c 'set -x; \
	time fastQDumpRunner.sh ${listOfSamples} ; \
	dataIsUpSMS.sh $phoneNumber $scriptName $fileOfSamples exit status: $? ' \
   2>&1 > $scriptName.${fileOfSamples}.${dateStamp}.out &



