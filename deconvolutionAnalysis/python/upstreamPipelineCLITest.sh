#
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
findModule="pipeline.dataFactory.test.exampleCreateSignatureGeneConfig"
estimatedScalingFactors="${rootDir}/pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/estimatedSizeFactors.csv"
outDir="${rootDir}/upstreamPipeline.out"
vargs="varg1 varg2 varg3"

rm -rf ${outDir}
#export PYTHONPATH="${PYTHONPATH}:${rootDir}"
python -m pipeline.upstreamPipeline \
       --colDataPath ${colData} \
       --countDataPath ${countData} \
       --deseqResultsDir ${deseqResultsDir} \
       --estimatedScalingFactors ${estimatedScalingFactors} \
       --findModule ${findModule} \
       --outDir ${outDir} \
       ${vargs}

