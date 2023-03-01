#!/bin/sh
# test to see if we can split a mixture matrix into parts and get same
# results as original mixture matrix
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
tmp=/scratch/aedavids/tmp
rootDir=/scratch/aedavids/GTEx_TCGA
bestRoot=geneSignatureProfiles/best
oneVsAll=GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25
#inputRootDir="${rootDir}/${bestRoot}"
bestSrc="${kl}/GTEx_TCGA/geneSignatureProfiles/best/${oneVsAll}"

mixtureMatrix=GTEx_TCGA_TrainGroupby_mixture.txt
'cp' -r "${bestSrc}/"ciberSort/${mixtureMatrix} ${tmp}

signatureMatrix=signatureGenes.tsv
'cp' -r "${bestSrc}/ciberSort/${signatureMatrix}" ${tmp}

#
# create 3 test mixture files
# cibersort mixtue format is genes x samples
#
cd ${tmp}

cut -f 1-100  ${mixtureMatrix} > mixtureMatrix100.tvs
cut -f 1-50   mixtureMatrix100.tvs > mixtureMatrixLeft50.tsv

cut -f 51-100 mixtureMatrix100.tvs > tmpMixtureMatrixRight50.tsv
cut -f 1 mixtureMatrix100.tvs > geneIds
paste geneIds tmpMixtureMatrixRight50.tsv > mixtureMatrixRight50.tsv

#
# configure docker
#

createJobId () {
    local mixtureFile=$1
    
    # dateStamp example: 2019-12-09-23.01.43-UTC
    local timeStamp=`date "+%Y-%m-%d-%H.%M.%S-%Z%n"`


    local jobId="${scriptName}-${mixtureFile}-${timeStamp}"

    # write return value to std out
    echo ${jobId}    
}

createOutputDir () {
    # # dateStamp example: 2019-12-09-23.01.43-UTC
    # local timeStamp=`date "+%Y-%m-%d-%H.%M.%S-%Z%n"`


    # local jobId="${scriptName}-${timeStamp}"
    local jobId=$1
    local outDir=/scratch/aedavids/cibersort.out/${jobId}
    mkdir -p $outDir

    # write return value to std out
    echo ${outDir}
}

# docker arguments
# -d  --detach Run container in background and print container ID
# -rm Automatically remove the container when it exits
# -e set environment variable
USER_ID=`id -u`
fmt="docker run \n\
    --detach \n\
    --rm \n\
    -e USERID=${USER_ID} \n\
    -v ${tmp}:/src/data \n\
    -v %s:/src/outdir \n\
    cibersortx/fractions \n\
    --username aedavids@ucsc.edu \n\
    --token 3f561ab6d4cf373d11f23d8e205b4b72 \n\
    --mixture %s \n\
    --sigmatrix ${signatureMatrix} \n\
    --perm 100 \n\
    --label %s \n\
    --QN FALSE \n\
    --verbose TRUE
"

#
# printf can be use like  sprintf() in 'c'
#


#mixture=mixtureMatrix100.tvs
#mixture=mixtureMatrixLeft50.tsv
mixture=mixtureMatrixRight50.tsv

jobId=`createJobId $mixture`
outDir=`createOutputDir $jobId`

cmd=$(printf "$fmt" "$outDir" "$mixture" "$jobId" )

# echo "********"
# printf $cmd
# echo "********"



printf "\n\n\n************ run docker\n"
echo $cmd > $outDir/${scriptName}.parameters.txt

scriptOut="$outDir/${scriptName}.meta.out"
echo "run on ${timeStamp}" > ${scriptOut}
#echo "input src: ${bestSrc}/ciberSort/*"  >> ${scriptOut}
echo $cmd >> ${scriptOut}
echo ""   >> ${scriptOut}


$cmd 
