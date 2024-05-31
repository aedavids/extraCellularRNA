#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 5/22/24
#
#
# comparing results from "best20GTEx_TCGA", "best25GTEx_TCGA", "best20FindAllDegree1_wl5", 
# suggest that finding degree 1 genes by looking at the largest number of deseq results
# rows should produce the best discriminators
#
# I do not expect we will add a lot of genes.
# the training time will be big
# next stage will remove genes from degree2 or greater intersections
# we can not remove these genes yet. I think it may cause some categories
# to drop out. There are 8 categories in best500 that do not have degreew
# 
# 
# runs analysis.OptimalSelectiveEnrichSignatureGeneConfig
#
# pipeline.sh is generic/reuseable.
# to reduce user and environment specifc coupling it takes a lot of arguments
# This script make calling the pipeline easier and makes the results reproducable
# we capture all the parameter values
#
# ref: extraCellularRNA/deconvolutionAnalysis/bin/test.slurm.sh
#

if [ $# -ne 2 ];
    then
        printf "ERROR  \n"
        printf "local install: cp ~/extraCellularRNA/deconvolutionAnalysis/bin/pipeline.sh .\n"
        printf "local install: cp ~/extraCellularRNA/deconvolutionAnalysis/bin/1vsAll-~gender_category/lumpBrain/${0} .\n"
        printf "usage: $0 ciberSortSecurityToken ciberSortUser\n"
        printf "follow 'Token and instruction access' @ https://cibersortx.stanford.edu/download.php"
        exit 1 # error
    fi

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

ciberSortSecurityToken=$1
ciberSortUser=$2

dataSet="GTEx_TCGA"

# rootDir="/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python"
rootDir="/private/groups/kimlab/${dataSet}"

# colData="${rootDir}/pipeline/dataFactory/test/data/testIntegration/colData.csv"
# colData="${rootDir}/groupbyGeneTrainingSets/${dataSet}_TrainColData.csv"
colData="${rootDir}/groupbyGeneTrainingSets/${dataSet}_TrainLumpBrain.colData.csv"

# countData="${rootDir}/pipeline/dataFactory/test/data/testIntegration/geneCounts.csv"
countData="${rootDir}/groupbyGeneTrainingSets/${dataSet}_TrainGroupby.csv"

# files in deseqResultsDir will be passed to 
# analysis.createOptimalSelectiveEnrichSignatureGeneConfig.findGenes()
# we want to enrich the results from the best500GTEx_TCGALumpBrain.sh stage
topN=500
designDir="1vsAll-~gender_category"
upstream="/private/groups/kimlab/aedavids/deconvolution/${designDir}/best${topN}${dataSet}LumpBrain/training"
upstreamOut="${upstream}/best${topN}${dataSet}LumpBrain.sh.out"
deseqResultsDir="${upstreamOut}/${dataSet}*"

# directories with results files we want to create historic gene set file from
historicGeneSetDirs="${deseqResultsDir}"

# search thes files for addtional degree1 genes
resultsDir="${rootDir}/1vsAll"

# findModule="analysis.createBest25CreateSignatureGeneConfig"
#findModule="analysis.createBestCreateSignatureGeneConfig"
#findModule="analysis.createBestRemoveHighDegreeSignatureGeneConfig"
#findModule="analysis.createEnrichSignatureGeneConfig"
findModule="analysis.createOptimalSelectiveEnrichSignatureGeneConfig"

# estimatedScalingFactors="${rootDir}/pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/estimatedSizeFactors.csv"
estimatedScalingFactors="${rootDir}/1vsAll/estimatedSizeFactors.csv"

outDir=`pwd`"/${0}.out"

# clean up any leftovers from old runs
'rm' -f serial*.log;

printf "deleted ${outDir}?\n"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) 'rm' -rf "${outDir}"; break;;
        No )  break;;
    esac
done

printf "deleted cromwell output?\n"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) 'rm' -rf cromwell-* ; break;;
        No )  break;;
    esac
done

'rm' -f import.zip
'rm' -f CIBERSORTxFractionsWorkflow.wdl
'rm' -f cromwellOptions.json 

mkdir -p "${outDir}"

printf "\n\n\n"

#
# prepare cromwell/wdl workflow input
# The pipeline will call cromwell and pass the cibersort input files
#
gitRoot="/private/home/aedavids/extraCellularRNA"
WDL_TOOLS="${gitRoot}/java/bin"


'cp' "${gitRoot}/deconvolutionAnalysis/bin/cromwellOptions.json.template" .

#
# copy required wdl files 
#
wdlRoot="${gitRoot}/terra/cibersortx/wdl"


rm -rf wdlTest
mkdir wdlTest
'cp' "${wdlRoot}/CIBERSORTxFractionsWorkflow.wdl" .
# -j just store file name ignore directory path
zip -j import.zip "${wdlRoot}/cibersortxFractionsTask.wdl" "${wdlRoot}/../wdlTest/partitionDataTask.wdl" "${wdlRoot}/../wdlTest/mergeTask.wdl"

wdlInputJSON="${outDir}/CIBERSORTxFractionsWorkflow.wdl.input.json"

#
# arguments to pass to OptimalSelectiveEnrichSignatureGeneConfig.__init__()
#

# BestN already found degree1 genes from idx 0. 
# we can improve performance by skipping these genes
startIdx=$topN

# search deseqResultsDir in range startIdx ... startIdx + windowLength
windowLength=500

# we want to find all degree1 genes
maxNumberOfGenes=10000

title=` printf "best%s_findAllDegree1_wl%s" $topN $windowLength` # do not uses spaces, title will be part of file paths

#
# we want to find all the degree1 genes in all categories
#
categories="Brain LIHC TGCT COAD Lung Esophagus_Gastroesophageal_Junction \
             Kidney_Cortex THCA Adrenal_Gland Breast_Mammary_Tissue \
            ACC KIRP LGG Uterus Artery_Coronary BRCA HNSC Nerve_Tibial \
             KICH  \
             Adipose_Subcutaneous  \
            CHOL Cells_Cultured_fibroblasts Cells_EBV-transformed_lymphocytes \
            PCPG Skin_Not_Sun_Exposed_Suprapubic Adipose_Visceral_Omentum Colon_Sigmoid \
             SKCM  DLBC Pituitary \
            Small_Intestine_Terminal_Ileum Vagina  Ovary MESO \
            Heart_Atrial_Appendage STAD Colon_Transverse Esophagus_Muscularis Whole_Blood \
            Artery_Aorta Pancreas THYM Cervix_Endocervix BLCA  \
             Skin_Sun_Exposed_Lower_leg LUSC LUAD PAAD CESC UCS GBM \
             Heart_Left_Ventricle  READ \
            Minor_Salivary_Gland Stomach PRAD Testis OV Esophagus_Mucosa Prostate \
            Artery_Tibial Muscle_Skeletal UCEC Bladder Spleen KIRC Thyroid SARC Liver ESCA UVM"

# example
# /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/
#       best20RemovedGTEx_TCGA/training/best20RemovedGTEx_TCGA.sh.out/
#       upsetPlot.out/best20_degreeThreshold_10.intersection.dict

#upsetOutBase="/private/groups/kimlab/aedavids/deconvolution/${designDir}/best${topN}${dataSet}/training"
upsetOutBase="/private/groups/kimlab/aedavids/deconvolution/${designDir}/best${topN}${dataSet}LumpBrain/training"

#upsetOut="${upsetOutBase}/best${topN}${dataSet}.sh.out/upsetPlot.out"
upsetOut="${upsetOutBase}/best${topN}${dataSet}LumpBrain.sh.out/upsetPlot.out"

#upStreamIntersectionDictionaryPath="${upsetOut}/best${topN}.intersection.dict"
upStreamIntersectionDictionaryPath="${upsetOut}/best${topN}.intersection.dict"

pythonSrcRoot="/private/home/aedavids/extraCellularRNA"

if [ -z ${PYTHONPATH+x} ];
    then
        #PYTHONPATH is unset or set to the empty string:
        export PYTHONPATH="${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
    else 
        export PYTHONPATH="${PYTHONPATH}:${pythonSrcRoot}/src:${pythonSrcRoot}/deconvolutionAnalysis/python"; 
    fi

printf "PYTHONPATH : $PYTHONPATH \n"

historicGeneSetPath="${outDir}/historicGeneSet.txt"
python -m pipeline.dataFactory.createHistoryGeneSet $historicGeneSetDirs > $historicGeneSetPath

vargs=" --design tilda_gender_category \
        --padjThreshold 0.001 \
        --lfcThreshold 2.0 \
        --dataSetName ${dataSet} \
        --windowLength ${windowLength} \
        --title ${title} \
        --localCacheRoot ${outDir} \
        --classes ${categories} \
        --historicGeneSetPath ${historicGeneSetPath} \
        --maxNumberOfGenes $maxNumberOfGenes \
        --resultsDir ${resultsDir} \
        --startIdx $startIdx \
        --upStreamIntersectionDictionaryPath ${upStreamIntersectionDictionaryPath} "


printf "\n\n\n SignatureGeneConfig vargs : $vargs \n !!!!!! \n"

logFile="${0}.log"
rm -f $logFile
setsid sh -c "set -x; pipeline.sh ${colData} \
                        ${countData} \
                        ${deseqResultsDir} \
                        ${findModule} \
                        ${estimatedScalingFactors} \
                        ${outDir} \
                        import.zip \
                        ${wdlInputJSON} \
                        ${WDL_TOOLS} \
                        ${gitRoot} \
                        ${ciberSortSecurityToken} \
                        ${ciberSortUser} \
                        ${vargs}" > $logFile 2>&1 & 

sleep 10
pstree $USER
 ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep $USER