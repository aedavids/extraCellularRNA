#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 12/28/2024
# 

scriptName=$0

if [ $# -ne 2 ];
    then
        scriptName=`basename $0`
        printf "ERROR  \n"
        printf "local install: cp ~/extraCellularRNA/deconvolutionAnalysis/tempus/bin/{${scriptName}} .\n"
        printf "usage: $0 ciberSortUser ciberSortSecurityToken \n"
        #printf "usage: tail -f ${0}.log \n"
        printf "follow 'Token and instruction access' @ https://cibersortx.stanford.edu/download.php \n"
        exit 1 # error
    fi

#
# Copy and past: Common paramters to change
# - topN
#

set -x
cibersortToken=$1
cibersortUser=$2



# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
#set -euxo pipefail
#set -eux pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

#
# select the biomarkers from the best10CuratedDegree1_ce467ff signatureGenes.tsv
#
sigDir='/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10CuratedDegree1_ce467ff/training/best10CuratedDegree1.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-10/ciberSortInput'
sigPath="${sigDir}/signatureGenes.tsv"

#  set -o pipefail and use head we get exit status 141 ????
rawHUGOSignatureGenes=`cut -f 1 ${sigPath} | sed -z 's/\n/ /g'`
#rawSignatureGenes=`cut -f 1 ${sigPath}| head `
#rawSignatureGenes=`cut -f 1 ${sigPath} | head`
#rawSignatureGenes=`cut -f 1 ${sigPath}| head | sed -e s/\n/ /g`
#rawSignatureGenes=`cut -f 1 ${sigPath} | head | xargs echo | sed -e 's/name//'`
#rawSignatureGenes="$(cut -f 1 ${sigPath} | head)"

# remove teh first entry, it is the header column name,
HUGOSignatureGenes=`echo ${rawHUGOSignatureGenes} | sed -z 's/name/ /g'`
printf "HUGOSignatureGenes: ${HUGOSignatureGenes} \n\n"

aedwip we need to mapp to V39


#
# set up the python env
#
printf "\n\n\n################ set up python environment \n"
# start conda env
condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`

# temporary turn +u off to prevent the following error
# miniconda3/envs/extraCellularRNA/etc/conda/activate.d/jgo_activate.sh: line 1: JGO_CACHE_DIR: unbound variable
#set +u # turn back off
source ${condaBase}/etc/profile.d/conda.sh
#set -u # turn back on
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

export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis";
printf "PYTHONPATH : $PYTHONPATH \n"

#
# create the cibersort mixture matrix
#

# cibersort input and output directory mount points
# You must use full paths. For unknown reasons cibersort raises an error if 
# the output directory in my home directory.
# the output and input directories can be the same
cibersortInputDir="${PWD}/${scriptName}.out/cibersortInputDir"
cibersortOutputDir="${PWD}/${scriptName}.out/cibersortOut" 
mkdir -p "${cibersortOutputDir}" "${cibersortOutputDir}"

normalizedCountFilePath="/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv"

printf "\n\n\n################### create the cibersort mixture matrix\n"
python -m tempus.createCibersortMixtureMatrix \
    --normalizedCountFilePath ${normalizedCountFilePath} \
    --genesOfInterest ${SignatureGenes} \
    --outDir ${cibersortInputDir}  

set -u # check for unset variables
#
# run the cibersortx docker container
#
printf "\n\n\n################### run the cibersortx docker container\n"
# ref : extraCellularRNA/terra/cibersortx/wdl/README.md
# ref : extraCellularRNA/terra/cibersortx/bin/run_cibersortx_fractions.sh
#aedwip cd ${cibersortInputDir}
mixtureMatrix=mixtureMatrix.tsv
signatureMatrix=signatureMatrix.tsv

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

exit 0


