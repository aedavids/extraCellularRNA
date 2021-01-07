#!/bin/bash

# aedavids@ucsc.edu


export scriptName=`basename $0`

if [ $# -ne 3 ]; then
    echo "error: usage $scriptName phoneNumber fileOfSamples salmonIndex"
    echo "example $scriptName 7508622639@txt.att.net SRR_Acc_List.txt.part.ad"
    echo "SRR_Acc_List.txt.part.ad is a list of samples"
    echo "fileOfSample creation example: 'split -l 100 SRR_Acc_List.txt SRR_Acc_List.txt.part.' "
    echo "example salmonIndex /public/groups/kimlab/indexes/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx"
    exit 1
fi

# if [ $# -lt 1 ]; then
#     printf "error: usage $scriptName salmonIndex dataDir outputDir\n"
#     printf "example $scriptName  \n\t /public/groups/kimlab/indexes/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx\  \n\t./data/PDAC/SRR10080544\ \n\t./data/PDAC/SRR10080544/salmon.out\n"
#     printf "assume dataDir has paired end data with file path format SRR10080544_pass_1.fastq.gz  SRR10080544_pass_2.fastq.gz \n"
#     printf "\n"
#     exit 1
# fi

# track arguments and environment setting to make more reproducable
echo cli: $0 $@
echo CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV
echo CONDA_PREFIX: $CONDA_PREFIX
echo `salmon --version`
echo 

# AEDWIP TODO improve packaging
# find a better place for scripts
which runSalmon.pancreas.plasma.ev.long.RNA.sh
exitStatus=$?
if [ ! $exitStatus -eq 0 ]; then
    echo "ERROR runSalmon.pancreas.plasma.ev.long.RNA.sh not found in PATH"
    echo " add ~/extraCellularRNA/bin"
    exit 1
fi

set -x # turn debug trace on
# set +x # turn debug trace off

# use export else variable will not be defined by inline script embedded in nohup
export phoneNumber=$1
export fileOfSamples=$2
export salmonIndex=$3
export listOfSamples=`cat ${fileOfSamples}`
export dateStamp=`dateStamp.sh`
export scriptLog="$scriptName.${fileOfSamples}.${dateStamp}.out"

# ugly hack
# make sure we are in the correct parent directory
export dataRoot=/public/groups/kimlab/pancreas.plasma.ev.long.RNA/
cd /public/groups/kimlab/pancreas.plasma.ev.long.RNA/

# create a session 
# this will ensure jobs continue to run even after we
# log out. It will also make easier to kill all the children
# if we need to restart for some reason
setsid sh -c 'set -x; \
        for sampleId in ${listOfSamples}; \
        do \
            dataDir=`ls -d ./data/*/${sampleId}`; \
            outputDir=${dataDir}/salmon.out; \
            log=${dataDir}/runSalmon.pancreas.plasma.ev.long.RNA.sh.${dateStamp}.out; \
            runSalmon.pancreas.plasma.ev.long.RNA.sh $salmonIndex $dataDir $outputDir 2>&1 > $log ;\
            printf "\n\n" ;\
        done ; \
        dataIsUpSMS.sh $phoneNumber $scriptName $fileOfSamples exit status: $? ' \
     > $scriptLog 2>&1 &


#


#             dataDir=`ls -d ./data/*/${sampleId}`; \
#             outputDir=${dataDir}/salmon.out; \
#             log=runSalmon.pancreas.plasma.ev.long.RNA.sh.${dateStamp}.out; \
# echo runSalmon.pancreas.plasma.ev.long.RNA.sh $salmonIndex $dataDir $outputDir 2>&1 > $log;\
