#!/bin/bash

#
# Andrew E. Davidson, aedavids@ucsc.edu
# 10/14/2021
#

set -x # turn trace debug log on
# set +x # turn trace debug log off

#R CMD BATCH hello.R

# CMD BATCH testDESeqScript.R --countMatrix mockCountMatrix.tsv --colData mockColData.tsv --design '~ + sampleType'

#R CMD testDESeqScript.R --countMatrix mockCountMatrix.tsv --colData mockColData.tsv --design '~ + sampleType' --referenceLevel zcontrol


# it does not matter if the model is '~ + foo' or '~ foo'
# R CMD testDESeqScript.R --countMatrix mockCountMatrix.tsv \
#   --colData mockColData.tsv \
#   --design '~ sampleType' \
#   --referenceLevel zcontrol \
#   --outFile "mockDESeq2_treatment_vs_control_results.tsv" \
#   --numCores 2

# R CMD DESeqScript.R \
#   --countMatrix masterCount.tsv/part-00000-d4d920b0-54d0-42bc-bc02-fa9cd1cde034-c000.csv \
#   --colData masterColData.tsv \
#   --design '~ treatment' \
#   --referenceLevel ctrl \
#   --outFile  masterDESeq2_kras_vs_control_results.tsv \
#   --numCores 2 \
#   --estimateSizeFactorsOutfile kras_vs_control_estimatedSizeFactors.tsv \
#   --oneVsAll

outDir="runner.sh.out"
mkdir -p $outDir
cd $outDir
R CMD ../DESeqScript.R \
  --countMatrix ../masterCount.tsv/part-00000-d4d920b0-54d0-42bc-bc02-fa9cd1cde034-c000.csv \
  --colData ../masterColData.tsv \
  --design '~ treatment' \
  --referenceLevel ctrl \
  --outFile  "masterDESeq2_kras_vs_control_results.tsv" \
  --numCores 2 \
  --estimateSizeFactorsOutfile "kras_vs_control_estimatedSizeFactors.tsv" \
  --oneVsAll


