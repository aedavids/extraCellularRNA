#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 1/1/24
# 
# ref : 
#   extraCellularRNA/terra/natureBioMedEng
#   extraCellularRNA/deconvolutionAnalysis/bin/pipeline.sh
#

#
# output env info to make debugging easier
#
pwd; 
hostname; 
date


# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

#
# clean up residuals from old runs
#
'rm' -rf cromwell-* ;

#
# set up the python env
# we do not use python, the conda environment has the correct version of java
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

#
# prepare cromwell/wdl workflow input
# The pipeline will call cromwell and pass the cibersort input files
#
gitRoot="/private/home/aedavids/extraCellularRNA"
WDL_TOOLS="${gitRoot}/java/bin"

wdlInputJSON="/private/groups/kimlab/aedavids/deconvolution/LUAD.vs.LUSC/LUAD.vs.LUSC.1vsAllTask.input.json"


# our conda environment should have the correct version of java
which java
java -version
printf "WDL_TOOLS : ${WDL_TOOLS}\n"


#
# start our wdl/docker runner
#      --imports "${wdlImportZip}" \
${gitRoot}/bin/runCromwell.sh \
     -Dconfig.file="${gitRoot}/terra/wdl/cromwellDebug.conf" \
     -jar "${WDL_TOOLS}/cromwell-85.jar" \
     run \
     --inputs "${wdlInputJSON}" \
     --options cromwellOptions.json \
     /private/home/aedavids/extraCellularRNA/terra/wdl/1vsAllTask.wdl

exitStatus=$?

echo ""
echo "runCromwell.sh exits status : $exitStatus"

exit $exitStatus

