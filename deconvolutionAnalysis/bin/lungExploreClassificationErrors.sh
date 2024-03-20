#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 12/08/23

numberOfArguments=4
if [ $# -lt $numberOfArguments ];
    then
        printf "ERROR missing command line arguments. expected $numberOfArguments recevied $# \n"
        echo "calculates falsePositives and falseNegatives for LUAD. Find shared genes between LUAD,LUSC,Lung,Whole_Blood"
        echo "usage: $0 intersectionDictionary results expected outDir"
        echo "ex: $0 analysis/test/data/intersection.dict analysis/test/data/results.tsv analysis/test/data/expectedFractions.tsv ${0}.out"
        exit 1 # error
    fi



# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

intersectionDictionary="$1"
results="$2"
expected="$3"
outDir="$4"

#
# set up the python env
#
printf "\n\n\nset up python environment \n"
# start conda env
condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`
source ${condaBase}/etc/profile.d/conda.sh
# set -x
conda activate extraCellularRNA

pythonSrcRoot="/private/home/aedavids/extraCellularRNA"

if [ -z ${PYTHONPATH+x} ];
    then
        #PYTHONPATH is unset or set to the empty string:
        export PYTHONPATH="${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
    else 
        export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
    fi

printf "PYTHONPATH : $PYTHONPATH \n"


python -m analysis.exploreClassificationErrors sg \
            --intersectionDictionary ${intersectionDictionary} \
            --category LUAD,LUSC,Lung,Whole_Blood \
            --outDir $outDir
  
for category in 'LUAD' 'LUSC' 'Lung' 'Whole_Blood' ;
do
    python -m analysis.exploreClassificationErrors fp \
            --category $category \
            --results $results  \
            --expected $expected \
            --outDir $outDir

    python -m analysis.exploreClassificationErrors fn \
            --category $category \
            --results $results  \
            --expected $expected \
            --outDir $outDir

done

echo "\n\n $0 finished"