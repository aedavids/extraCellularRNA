#!/bin/sh
# run cibersort docker on GTEX_TCGA train data
#
# Andrew E. Davidson
# aedavids@ucsc.edu
# 


scriptName=`basename $0`

# exit when any command fails
set -e

set -x   # turn debug on
# set +x # turn debug off


#
# copy ciber sort input from /private (nfs) to /scratch (local hard disk)
#
rootDir=/scratch/aedavids/GTEx_TCGA
bestRoot=geneSignatureProfiles/best
oneVsAll=GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25
inputRootDir="${rootDir}/${bestRoot}"

# cp and docker do not work with file names containing non alpha numeric chars
fixAlphaNumeric="${oneVsAll}/ciberSort"
inputDir="${inputRootDir}/${fixAlphaNumeric}"
mkdir -p ${inputDir}

#bestSrc="${kl}/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25"
bestSrc="${kl}/GTEx_TCGA/geneSignatureProfiles/best/${oneVsAll}"

#
# The non-alpha numeric chars in the file path cause problems
# docker can not mount volumns. cp does not work with out stranging quoting
# symbolic link syntax ls -s target directory
#
hack=${inputRootDir}/tmp
if [ ! -L "$hack" ]; then
    # if link exits ln will generate error causing script to end
    ln -s ${inputDir} ${hack}
fi

inputDir=${hack}
cp -r "${bestSrc}/"ciberSort/* ${inputDir}

printf "\n\n\n********** create output dir\n"

#
# create output directory on /scratch
#

# dateStamp example: 2019-12-09-23.01.43-UTC
timeStamp=`date "+%Y-%m-%d-%H.%M.%S-%Z%n"`


jobId="GTEx_TCGA_TrainGroupby_mixture-${timeStamp}"
outDir=/scratch/aedavids/cibersort.out/${jobId}
mkdir -p $outDir

#
# run docker
# create cmd string we can use to capture the optional parameters
#
printf "\n\n\n************ configure docker\n"
cd ${inputDir}
mixtureMatrix=GTEx_TCGA_TrainGroupby_mixture.txt
signatureMatrix=signatureGenes.tsv


# docker arguments
# -d  --detach Run container in background and print container ID
# -rm Automatically remove the container when it exits
# -e set environment variable

USER_ID=`id -u`
cmd="docker run \
    --detach \
    --rm \
    -e USERID=${USER_ID} \
    -v ${inputDir}:/src/data \
    -v ${outDir}:/src/outdir cibersortx/fractions \
    --username aedavids@ucsc.edu \
    --token 3f561ab6d4cf373d11f23d8e205b4b72 \
    --mixture ${mixtureMatrix}\
    --sigmatrix ${signatureMatrix}\
    --perm 100 \
    --label $jobId \
    --QN FALSE \
    --verbose TRUE
"

printf "\n\n\n************ run docker\n"
echo $cmd > $outDir/${scriptName}.parameters.txt

scriptOut="$outDir/${scriptName}.meta.out"
echo "run on ${timeStamp}" > ${scriptOut}
echo "input src: ${bestSrc}/ciberSort/*"  >> ${scriptOut}
echo $cmd >> ${scriptOut}
echo ""   >> ${scriptOut}

echo ""
$cmd 
