#!/bin/bash
# aedavids@ucsc.edu
# https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/


# run fastqc on a couple of file to get rough idea about quality


echo cli: $0 $@
echo CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV
echo CONDA_PREFIX: $CONDA_PREFIX
echo `fastqc --version`

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
Healthy_Control="Healthy_Control/Sample_N1/Sample_N1.fastq.zip Healthy_Control/Sample_N2/Sample_N2.fastq.zip
"
Colon_Carcinoma="Colon_Carcinoma/Sample_1S1/Sample_1S1.fastq.zip"
Pancreatic_Carcinoma="Pancreatic_Carcinoma/Sample_Pan01/Sample_Pan01.fastq.zip"
Prostate_Carcinoma="Prostate_Carcinoma/Sample_PC1/Sample_PC1.fastq.zip"

for i in $Healthy_Control $Colon_Carcinoma $Pancreatic_Carcinoma $Prostate_Carcinoma;
do
    echo "" # add white space to make log output easier to read
    prefix=`echo $i | cut -d / -f 1,2` # Colon_Carcinoma/Sample_1S1
    #echo $prefix

    outdir="${dataRoot}/${prefix}/fastqc.out"
    mkdir -p $outdir
    
    originalDataFile="${dataRoot}/${prefix}/*.zip"
    trimmedDataFile="${dataRoot}/${prefix}/trimmomatic.out/*.gz"

    set -x # turn debug trace on

    # hack to get fastqc to read a zip file
    # https://github.com/s-andrews/FastQC/issues/14
    # name=`echo $prefix | tr "/" ":"`
    # zcat $originalDataFile | fastqc --outdir $outdir --thread 6 "stdin:${name}"
    fastqc --outdir $outdir --thread 6 ${trimmedDataFile}

    set +x # turn debug trace off 
done
