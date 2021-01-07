#!/bin/bash

# aedavids@ucsc.edu
# hacky little script to move /public/groups/kimlab/pancreas.plasma.ev.long.RNA/fastq
# into sub directories

#set -x # turn debug on
# set +x # turn debug off

echo running $0 $@
echo `fastq-dump --version`

rootDir=/public/groups/kimlab/pancreas.plasma.ev.long.RNA
metaDataFile=${rootDir}/SraRunTable.csv
fastQDir=${rootDir}/fastq

# sed 1d SraRunTable.csv $ remove the header line from the csv
# cut 1 = sample Id , 19 = disease_state
metaData=`sed 1d ${metaDataFile} | cut -d , -f 1,19`

for line in $metaData
do
    sampleId=`echo $line | cut -d, -f 1`
    diseaseState=`echo $line | cut -d, -f 2`
    targetDir="${rootDir}/data/${diseaseState}/${sampleId}/"
    echo # add white space for readablity 
    set -x # turn debug on
    # set +x # turn debug off

    mkdir -p ${targetDir}
    mv $fastQDir/${sampleId}* ${targetDir}

    #ls $fastQDir/${sampleId}*
    #set -x # turn debug on
    set +x # turn debug off

    # add white space for readablity
    echo


done

