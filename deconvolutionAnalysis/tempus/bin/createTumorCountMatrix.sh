#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 1/14/2025
# 
# reads tempus/illumina/20241107/create/results/data/raw_counts.csv
# these are the group by gene counts from the salmon quant.sf files
# each row is a sample and the columns are the gene ids
# 
# output
# for each tumor type
#   1. create createTumorCountMatrix.sh.output
#   2. create tumor id RawGroupByGenesCounts.csv. each row is a sample and the columns are the gene ids
#   3. create tumor id RawDESeq2FmtCounts.csv. This is the transpose of the RawGroupByGenesCounts.csv
#   they are in the format DESeq expects 
#   4. create tumor id RawDESeq2FmtControlUDCounts.csv.csv. This is the transpose of the RawGroupByGenesCounts.csv
#      with only the Control and UD samples in the format DESeq@ expects



set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

#
# output env info to make debugging easier
#
scriptName=`basename $0`
printf "\n ######### script name: %s BEGIN \n" $scriptName
pwd; 
hostname; 
date

outDir="`pwd`/${scriptName}.output"
mkdir -p "${outDir}"

dataDir=/private/groups/kimlab/data/tempus/illumina/20241107/create/results/data

#
# set up the python env
#
printf "\n\n\n################ set up python environment \n"
# start conda env
condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`
source ${condaBase}/etc/profile.d/conda.sh
# set -x
conda activate extraCellularRNA

pythonSrcRoot="/private/home/aedavids/extraCellularRNA"

tempusPythonPath="${pythonSrcRoot}/deconvolutionAnalysis/python:${pythonSrcRoot}/deconvolutionAnalysis/python/tempus"
if [ -z ${PYTHONPATH+x} ];
    then
        #PYTHONPATH is unset or set to the empty string:
        export PYTHONPATH="${pythonSrcRoot}/src:${tempusPythonPath}"; 
    else 
        export PYTHONPATH="${PYTHONPATH}:${tempusPythonPath}"; 
    fi

printf "\n ########### PYTHONPATH : $PYTHONPATH \n"


#
# raw_counts.csv values are the group by gene counts from the salmon quant.sf files
# the last column is the sample id. all the other columns names are gene ids.
#
# $ head -n 1 raw_counts.csv | comma2newLine | wc -l
# 76540
#
# $ head raw_counts.csv | cut -d , -f 1,2,3,4,76539,76540
# (A)n,(AAA)n,(AAAAAAC)n,(AAAAAAG)n,Zaphod3,gene_id
# 1,0,0,0,9,SLDK3_T1_100_S1_L007
# 2,0,0,0,0,SLDK3_T1_100K_S2_L007
# 1,0,0,0,21,SLDK3_T1_10K_S3_L007
# 0,0,0,0,0,SLDK3_T1_1K_S4_L007
# 1,0,0,0,0,SLDK3_T1_1M_S5_L007
# 0,0,0,0,0,SLDK3_T1_Control_S6_L007
# 5,0,0,0,0,SLDK3_T1_UD_S7_L007
# 0,0,0,0,1,SLDK3_T2_100_S8_L007
# 0,0,0,0,6,SLDK3_T2_100K_S9_L007

# tumorId=T1
# tumorToken="_${tumorId}_"
# grep $tumorToken "${dataDir}/raw_counts.csv" > "${outDir}/${tumorId}GroupByGenesCounts.csv"

# ?? T4 is missing
validTokens=("T1" "T2" "T3" "T5" "T6")

for tumorId in "${validTokens[@]}"; 
do
    # printf "tumorId: %s\n" $tumorId
    tumorToken="_${tumorId}_"

    # copy first line of raw file into output file
    # it has the gene ids
    outRawFile="${outDir}/${tumorId}RawGroupByGenesCounts.csv"
    head -n 1 "${dataDir}/raw_counts.csv" > "${outRawFile}"

    # copy all the lines with the tumor token into the output file
    grep $tumorToken "${dataDir}/raw_counts.csv" >> "${outRawFile}"
    # printf "_${tumorId}_ exit code $? \n"

    outAllFile="${outDir}/${tumorId}RawDESeq2FmtCounts.csv"
    controlUDFile="${outDir}/${tumorId}RawDESeq2FmtControlUDCounts.csv"
    
    #run python to transpose
    python -m tempus.transposeCounts \
        --tumorCountMatrixFilePath "${outRawFile}" \
        --outAllPath "${outAllFile}" \
        --outCntrlUDPath "${controlUDFile}"

done

printf "\n ######### script name: %s END \n" $scriptName

