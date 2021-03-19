#!/bin/bash

# aedavids@ucsc.edu

#
# select the instDir
# select the batchList
#


set -x # turn trace debug log on

#instDir=donorMetaInstr
#instDir=fastQBatch
#instDir=bioSampleMetaInstr
instDir=exprMetaInstr

rootDir=/public/groups/kimlab/plasma.ex.RNA.Patel
instrRoot="${rootDir}/downloadInstructions/${instDir}"
dataRoot="${rootDir}/data"

filePrefix="${instDir}.a"

# run 8 downloads at a time. each instruction files downloads about 10 files
batchList1="a b c d e f g" # donorMetaInstr fastQBatch exprMetaInstr
batchList2="h i j k l m n" # donorMetaInstr fastQBatch exprMetaInstr
batchList3="p q r s t"     # donorMetaInstr fastQBatch exprMetaInstr bioSampleMetaInstr

batchList=$batchList1
#batchList=$batchList2
#batchList=$batchList3
for i in $batchList
do
	downLoadInstrFile="${instrRoot}/${filePrefix}$i"
	batchDownLoadExRNA.orgData.sh ${downLoadInstrFile} ${dataRoot} 6508622639@txt.att.net
done
