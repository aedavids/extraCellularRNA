#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 12/28/2024
# 

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
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

# where if we use set -o pipe fail we get exit status 141 ????
rawSignatureGenes=`cut -f 1 ${sigPath} | sed -z 's/\n/ /g'`
#rawSignatureGenes=`cut -f 1 ${sigPath}| head `
#rawSignatureGenes=`cut -f 1 ${sigPath} | head`
#rawSignatureGenes=`cut -f 1 ${sigPath}| head | sed -e s/\n/ /g`
#rawSignatureGenes=`cut -f 1 ${sigPath} | head | xargs echo | sed -e 's/name//'`
#rawSignatureGenes="$(cut -f 1 ${sigPath} | head)"

# remove teh first entry, it is the header column name,
SignatureGenes=`echo ${rawSignatureGenes} | sed -z 's/name/ /g'`
printf "SignatureGenes: ${SignatureGenes} \n"


# #
# # set up the python env
# #
# printf "\n\n\n################ set up python environment \n"
# # start conda env
# condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`
# source ${condaBase}/etc/profile.d/conda.sh
# # set -x
# conda activate extraCellularRNA

# pythonSrcRoot="/private/home/aedavids/extraCellularRNA"

# if [ -z ${PYTHONPATH+x} ];
#     then
#         #PYTHONPATH is unset or set to the empty string:
#         export PYTHONPATH="${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
#     else 
#         export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
#     fi

# printf "PYTHONPATH : $PYTHONPATH \n"

# #
# # create the cibersort mixture matrix
# #
# printf "\n\n\n################### create the cibersort mixture matrix\n"
# python -m tempus.createCibersortMixtureMatrix \
#     --normalizedCountFilePath ${normalizedCountFilePath} \
#     --genesOfInterest ${genesOfInterest} \
#     --outDir ${cibersortInputDir}  

# exit 0