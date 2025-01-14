#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 1/14/2025
# 
# create a count matrix from

set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

#
# output env info to make debugging easier
#
scriptName=`basename $0`
pwd; 
hostname; 
date

# the salmon output files.
dataDir=/private/groups/kimlab/data/tempus/illumina/20241107/create/quant

# valid tokens are T1, T2, T3, T4, T5, T6
AEDWIP_TUMOR_TOKEN=T1
tumorToken=$AEDWIP_TUMOR_TOKEN

tumorSamples=`ls ${dataDir} | grep ${tumorToken}`


for i in $tumorSamples;
do
    count=0;
    for j in $tumorSamples;
    do
        #c=$(($a + $b))
        # x=1
        count=$(($count + 1))
        # count=$c
        printf "count : $count\n"
        numReadsColPosition=5
        salmonResults="${dataDir}/${j}/quant.sf"
        if [ $count -eq 1 ]; then
            echo "This is the first iteration."
        fi
        head $salmonResults | cut -f $numReadsColPosition
    done

    printf "end i : $i *************\n"
done