#!/bin/bash
# aedavids@ucsc.edu
# https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/


#
# hack to see if we need to remove adapters or not
# assume you are  the /public/groups/kimlab/pancreas.plasma.ev.long/fastq


scriptName=`basename $0`
# if [ $# -lt 1 ]; then
#     echo "error: usage $scriptName  list of SRA tool sample ids"
#     echo "example $scriptName  SRR10080507*"
#     exit 1
# fi


echo cli: $0 $@
echo CONDA_DEFAULT_ENV: $CONDA_DEFAULT_ENV
echo CONDA_PREFIX: $CONDA_PREFIX
echo `fastqc --version`

set -x # turn debug trace on
# set +x # turn debug trace off

hack=/public/home/aedavids/extraCellularRNA/bin
now=`$hack/dateStamp.sh`
outdir=fastqc.${now}.out
mkdir $outdir
adapterDir=/public/groups/kimlab/genomes.annotations/adapters
#adapter=$adapterDir/TruSeq-HT.fastqc.txt
adapter=$adapterDir/TruSeq-HT.fastqc.tab.txt
#fastqc --outdir $outdir --adapters $adapter --thread 6 SRR9624972* SRR9624860* SRR9624849*

fastqc --outdir $outdir --adapters $adapter --thread 6 SRR9625015* SRR9625003*
