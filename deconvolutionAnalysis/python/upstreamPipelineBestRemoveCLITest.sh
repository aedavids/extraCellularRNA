#
# test analysis/bestRemoveHighDegreeSignatureGeneConfig.py 
# Andrew E. Davidson, aedavids@ucsc.edu
# 11/14/23
#

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

pwd

rootDir=`dirname $0`

colData="${rootDir}/pipeline/dataFactory/test/data/testIntegration/colData.csv"
countData="${rootDir}/pipeline/dataFactory/test/data/testIntegration/geneCounts.csv"
deseqResultsDir="${rootDir}/pipeline/dataFactory/test/data/testSignatureGenes/1vsAll" 
findModule="analysis.createBestRemoveHighDegreeSignatureGeneConfig"
estimatedScalingFactors="${rootDir}/pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/estimatedSizeFactors.csv"
outDir="${rootDir}/${0}.out"

# arguments to pass to BestRemoveHighDegreeSignatureGeneConfig.__init__()
topN=20
title=` printf "best%s" $topN` # do not uses spaces, title will be part of file paths
vargs=" --design tilda_gender_category \
        --padjThreshold 0.001 \
        --lfcThreshold 2.0 \
        --dataSetName GTEx_TCGA \
        --number  $topN \
        --title ${title} \
        --localCacheRoot ${outDir} \
        --intersectionDict analysis/test/data/intersection.dict \
        --degreeThreshold 1"

rm -rf ${outDir}
#export PYTHONPATH="${PYTHONPATH}:${rootDir}"
python -m pipeline.upstreamPipeline \
       --colDataPath ${colData} \
       --countDataPath ${countData} \
       --deseqResultsDir ${deseqResultsDir} \
       --estimatedScalingFactors ${estimatedScalingFactors} \
       --findModule ${findModule} \
       --outDir ${outDir} \
       --vargs ${vargs}
