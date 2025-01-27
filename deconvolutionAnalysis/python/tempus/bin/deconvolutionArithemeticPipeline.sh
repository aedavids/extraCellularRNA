#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 12/25/2024
# 

#
# output env info to make debugging easier
#
scriptName=`basename $0`
pwd; 
hostname; 
date

# print all the cli argumnents
printf "################ ${scriptName} : command line arguments \n"
for var in "$@"
do
    printf "argument :$var \n"
done

# parse the arguments
numberOfArguments=8 

if [ $# -ne $numberOfArguments ];
    then
        printf "ERROR ${scriptName} missing command line arguments. expected $numberOfArguments recevied $# \n"
        printf "$@ \n"
        exit 1 # error
    fi

set -x

# cibersortx arguments
cibersortUser=$1
cibersortToken=$2

#tempus.findArithmeticBiomarkers arguments
topN=$3
normalizedCountsPath="${4}"
metaDataFilePath="${5}"
tumorId=$6

# cibersort mount points must be full paths
cibersortInputDir="${7}"
cibersortOutputDir="${8}"
mkdir -p "${cibersortInputDir}" "${cibersortOutputDir}"

set -x turn debug trace on. output goes to stderr, normal output goes to stdout

set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options


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

if [ -z ${PYTHONPATH+x} ];
    then
        #PYTHONPATH is unset or set to the empty string:
        export PYTHONPATH="${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
    else 
        export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
    fi

printf "PYTHONPATH : $PYTHONPATH \n"


printf "\n\n\n################### findArithmeticBiomarkers\n"
python -m tempus.findArithmeticBiomarkers \
    --normalizedCountsPath ${normalizedCountsPath} \
    --outDir ${cibersortInputDir}  \
    --topN ${topN} \
    --tumorId ${tumorId}


#
# run the cibersortx docker container
#
printf "\n\n\n################### run the cibersortx docker container\n"
# ref : extraCellularRNA/terra/cibersortx/wdl/README.md
# ref : extraCellularRNA/terra/cibersortx/bin/run_cibersortx_fractions.sh
#aedwip cd ${cibersortInputDir}
mixtureMatrix=mixtureMatrix_${tumorId}.tsv
signatureMatrix=signatureMatrix_${tumorId}.tsv

# dateStamp example: 2019-12-09-23.01.43-UTC
timeStamp=`date "+%Y-%m-%d-%H.%M.%S-%Z%n"`
jobId="${scriptName}-${timeStamp}"

# docker arguments
# -d  --detach Run container in background and print container ID
# -rm Automatically remove the container when it exits
# -e set environment variable

# mount the directory with the signature and matrix files as /src/data
# use full path 
# weird for unknown reasons cibersort raises an error if the output directory
# in my home directory.
# cibersortOutputDir="${PWD}/cibersortOutputDir"
# cibersortOutputDir=/scratch/aedavids/cibersortOut 
# mkdir -p $cibersortOutputDir

img="cibersortx/fractions"

printf "\n\n\n************ AEDWIP do not run docker we do not have a valid token\n"
#
# you can debug docker token problems by removing
# the --detach flag, --rm and adding --interactive --tty
# this will cause container error message to be written to the terminal
#
#    --detach \
#     --rm \
USER_ID=`id -u`
cmd="docker run \
   --detach \
    --rm \
    -e USERID=${USER_ID} \
    -v ${cibersortInputDir}:/src/data \
    -v ${cibersortOutputDir}:/src/outdir \
    ${img} \
    --username ${cibersortUser} \
    --token ${cibersortToken} \
    --mixture ${mixtureMatrix}\
    --sigmatrix ${signatureMatrix}\
    --perm 100 \
    --label $jobId \
    --QN FALSE \
    --verbose TRUE
"

printf "\n\n\n************ run docker\n"
echo "${cmd}" > "${cibersortInputDir}/${scriptName}.parameters.txt"

# scriptOut="$cibersortInputDir/${scriptName}.meta.out"
# echo "run on ${timeStamp}" > ${scriptOut}
# echo "input src: ${bestSrc}/ciberSort/*"  >> ${scriptOut}
# echo $cmd >> ${scriptOut}
# echo ""   >> ${scriptOut}

echo ""
$cmd 
exitStatus=$?
if [ $exitStatus -ne 0 ]; then
    printf "docker run failed with exit status : ${exitStatus}\n"
    exit $exitStatus
fi
