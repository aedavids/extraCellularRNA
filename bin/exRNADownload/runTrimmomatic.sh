#!/bin/bash
# aedavids@ucsc.edu
# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf


# run trimmomatic to see if it resolve quality issues identified
# using fastqc on a couple of file


echo cli: $0 $@
echo CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV
echo CONDA_PREFIX: $CONDA_PREFIX
echo trimmomatic version `trimmomatic -version`

# set -x # turn debug trace on
# set +x # turn debug trace off 

root="/public/groups/kimlab/plasma.ex.RNA.Patel"
if [ `pwd` != $root ]; then
    echo ERROR script hack you must be in  $root
    exit 1
fi

dataRoot="${root}/data"

#
# data files
#
Healthy_Control="Healthy_Control/Sample_N1/Sample_N1.fastq.zip 
                 Healthy_Control/Sample_N2/Sample_N2.fastq.zip"

Colon_Carcinoma="Colon_Carcinoma/Sample_1S1/Sample_1S1.fastq.zip"

Pancreatic_Carcinoma="Pancreatic_Carcinoma/Sample_Pan01/Sample_Pan01.fastq.zip"

Prostate_Carcinoma="Prostate_Carcinoma/Sample_PC1/Sample_PC1.fastq.zip"

#
# configure ILLUMINACLIP parameters
#
adapterRoot="/public/home/aedavids/miniconda3/envs/extraCellularRNA/share/trimmomatic-0.39-1/adapters"
#adapter="${adapterRoot}/TruSeq2-SE.fa"
adapter="${adapterRoot}/TruSeq3-SE.fa" # illumina Universal Adapter


# Roman's paired parameters
# :1:30:10:4:true 
#
# SE examplehttps://angus.readthedocs.io/en/2019/quality-and-trimming.html
# http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
# "short sections of each adapter (maximum 16 bp) are tested in each possible position
#  within the reads. If this short alignment, known as the 'seed' is a perfect or
#  sufficiently close match, determined by the seedMismatch parameter (see below), 
#  the entire alignment between the read and adapter is scored"
seedMismatches="1"

# palindromeClipThreshold: specifies how accurate the match between the two 
# 'adapter ligated' reads must be for PE palindrome read alignment
# palindromeClipThreshold="30"
# use zero because we do not have forward and reverse reads
palindromeClipThreshold="0"

# simpleClipThreshold: specifies how accurate the match between any adapter etc.
# sequence must be against a read
simpleClipThreshold="10"


# optional which affect palindrome mode only
# <minAdapterLength>:<keepBothReads>

illuminaClipParams="${adapter}:${seedMismatches}:${palindromeClipThreshold}:${simpleClipThreshold}"

for i in $Healthy_Control $Colon_Carcinoma $Pancreatic_Carcinoma $Prostate_Carcinoma;
do
    prefix=`echo $i | cut -d / -f 1,2` # Colon_Carcinoma/Sample_1S1
    #echo $prefix

    outDir="${dataRoot}/${prefix}/trimmomatic.out"
    mkdir -p $outDir

    echo ""
    set -x # turn debug trace on    

    # convert zip to gz
    zipDataFile="${dataRoot}/${prefix}/*.zip"
    gzInputDataFile=`echo $zipDataFile | sed s/fastq.zip/fastq.gz/`
    if [ ! -f ${gzInputDataFile} ]; then
        zcat ${zipDataFile} | gzip -c > ${gzInputDataFile}
    fi

    # output gz file
    outFileName=`basename $gzInputDataFile  | sed s/fastq.zip/fastq.gz/`
    gzOutFile="${outDir}/${outFileName}"

    logFile="${outDir}/trimmomatic.log"

    summaryLog="${outDir}/summary.log"

    # processing of different steps occur in ocer in which the steps
    # are specified on cli
    trimmomatic SE -threads 16 \
        -trimlog $logFile \
        -summary $summaryLog \
        $gzInputDataFile $gzOutFile \
        ILLUMINACLIP:${illuminaClipParams}


    # clean up
    # TODO AEDWIP while we are debugging leave the gz version of the original
    # will make debugging faster
    #'rm' ${gzInputDataFile}

    set +x # turn debug trace off 


done

