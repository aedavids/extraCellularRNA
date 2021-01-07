#!/bin/bash

# aedavids@ucsc.edu
# ref: /public/home/rreggiar/projects/aale.kras/scripts/salmonRun.sh

scriptName=`basename $0`
if [ $# -lt 1 ]; then
    printf "error: usage $scriptName salmonIndex dataDir outputDir\n"
    printf "example $scriptName  \n\t /public/groups/kimlab/indexes/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx\  \n\t./data/PDAC/SRR10080544\ \n\t./data/PDAC/SRR10080544/salmon.out\n"
    printf "assume dataDir has paired end data with file path format SRR10080544_pass_1.fastq.gz  SRR10080544_pass_2.fastq.gz \n"
    printf "\n"
    exit 1
fi

# track arguments and environment setting to make more reproducable
echo cli: $0 $@
echo CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV
echo CONDA_PREFIX: $CONDA_PREFIX
echo `salmon --version`
echo 
set -x # turn debug trace on
# set +x # turn debug trace off

salmonIndex="$1"
dataDir="$2"
outputDir="$3"


if [[ ! -f "$outputDir"/quant.sf ]]; then

    mkdir -p "$outputDir"
    
    #set +x # turn debug trace off

    sampleId=`basename ${dataDir}`
    trim_fwd="${dataDir}/${sampleId}_pass_1.fastq.gz"
    trim_rev="${dataDir}/${sampleId}_pass_2.fastq.gz"
    #set -x # turn debug trace on

    # https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-mapping-based-mode
    # --libType A : automatically infer the library type
    # --gcBias: model will attempt to correct for biases in how likely a
    #      sequence is to be observed based on its internal GC content
    # --seqBias:  will attempt to correct for random hexamer priming bias, 
    #          which results in the preferential sequencing of fragments
    #          starting with certain nucleotide motifs.

    # AEDWIP  --recoverOrphans : only be used in conjunction with selective alignment)

    salmon quant \
        -i "$salmonIndex" \
        --libType A \
        -1 "$trim_fwd" \
        -2 "$trim_rev" \
        -p 8 \
        --recoverOrphans \
        --validateMappings \
        --gcBias \
        --seqBias \
        --rangeFactorizationBins 4 \
        --output "$outputDir" 

    exitStatus=$?
    #set -x # turn debug on
    #set +x # turn debug off

    if [ $? -ne 0 ]; then

        echo ERROR salmon "$dataDir" returned exit status "$exitStatus"
        continue

    fi

fi    
#}



