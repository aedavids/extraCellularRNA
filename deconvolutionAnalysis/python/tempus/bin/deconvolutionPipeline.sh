#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 12/25/2024
# 

# TODO AEDWIP parse command line arguments
scriptName=`basename $0`

# tempus.findBiomarkers arguments
lfcThreshold=2.0 
topN=10
padjThreshold=0.001
deseq2ResultsFilePath="/private/groups/kimlab/data/tempus/illumina/20241107/create/results/data/Control_vs_UD_deseq_results.csv"

# cibersort mount points must be full paths
#cibersortInputDir="./tmp"
cibersortInputDir="${PWD}/cibersortInputDir"

# tempus.createCibersortSignatureMatrix arguments
#TODO add useMedian argument
normalizedCountFilePath="/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv"
metaDataFilePath="/private/groups/kimlab/data/tempus/illumina/20241107/raw/metaDataWithHeader.csv"
categoriesOfInterest="UD Control"

# cibersortx arguments
cibersortUser="aedavids@ucsc.edu"
cibersortToken="aac8366d23af037ef2423d0dbe8fffd8"

set -x turn debug trace on. output goes to stderr, normal output goes to stdout

set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options


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


#
# find biomarkers
# 
python -m tempus.findBiomarkers  \
    --lfcThreshold ${lfcThreshold} \
    --number ${topN} \
    --padjThreshold ${padjThreshold} \
    --deseq2ResultsFilePath ${deseq2ResultsFilePath} \
    --outDir ${cibersortInputDir} 

#
# get the top N biomarkers
# use cut to get the gene_id column
# grep -v remove the column header
# xargs echo -n remove the newline
#
genesOfInterest=`cat ${cibersortInputDir}/biomarkerDESeq2Results.csv | cut -d , -f 7 | grep -v gene_id | xargs echo -n`
printf "genesOfInterest : ${genesOfInterest}\n"

#
# create the signature matrix for cibersortx
#
#TODO add useMedian argument
printf "\n\n\n################### create the signature matrix for cibersortx\n"
python -m tempus.createCibersortSignatureMatrix \
    --normalizedCountFilePath ${normalizedCountFilePath} \
    --genesOfInterest ${genesOfInterest} \
    --metaDataFilePath ${metaDataFilePath} \
    --outDir ${cibersortInputDir}  \
    --categoriesOfInterest ${categoriesOfInterest}

#
# create the cibersort mixture matrix
#
printf "\n\n\n################### create the cibersort mixture matrix\n"
python -m tempus.createCibersortMixtureMatrix \
    --normalizedCountFilePath ${normalizedCountFilePath} \
    --genesOfInterest ${genesOfInterest} \
    --outDir ${cibersortInputDir}  


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
cibersortOutputDir="${PWD}/cibersortOutputDir"
cibersortOutputDir=/scratch/aedavids/cibersortOut 
mkdir -p $cibersortOutputDir

img="cibersortx/fractions"

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
