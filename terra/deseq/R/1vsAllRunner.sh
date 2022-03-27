#!/bin/bash

#
# Andrew E. Davidson, aedavids@ucsc.edu
# 1/25/2022
#

scriptName=$0
set -x # turn trace debug log on
# set +x # turn trace debug log off

# it does not matter if the model is '~ + foo' or '~ foo'


RSCRIPT="DESeqScript.R"
# outdir is docker mount point. assumes we are running container on the machine mustard
rootDir="/scratch/aedavids/GTExData"
outdir="${rootDir}/1vsAllRunner.sh.out"
mkdir -p "${outdir}"
#refLevel="Lung"
#refLevel="Pancreas"
refLevel="Thyroid"
dataSet="Validate"
outPrefix="${outdir}/${dataSet}_${refLevel}_vs_all"
#cm="${rootDir}/cleanCountsGroupedByGene/part-00000-1de2d5ed-05bd-4bbc-8cbf-ac2a3f8b98b2-c000.csv "
cm="${rootDir}/data/deseq/GTEx${dataSet}GroupByGenesCountMatrix.csv"
R CMD $RSCRIPT \
  --countMatrix $cm \
  --colData "${rootDir}/data/deseq/GTEx${dataSet}ColData.csv" \
  --design '~ sex + tissue_id' \
  --referenceLevel ${refLevel} \
  --outFile  "${outPrefix}_results.csv" \
  --numCores 6 \
  --estimateSizeFactorsOutfile "${outPrefix}_estimatedSizeFactors.csv" \
  --oneVsAll \
  --isCSV 2>&1 > ${scriptName}.out



echo " "
echo "*********** check ${RSCRIPT}.out for [R] stdout"
