#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 11/15/23
# 
# to use slurm see extraCellularRNA/deconvolutionAnalysis/bin/test.slurm.sh
#
# ref:
#   https://giwiki.gi.ucsc.edu/index.php/Overview_of_using_Slurm
#   extraCellularRNA/terra/cibersortx/wdl/README.md 
#   extraCellularRNA/terra/cibersortx/wdl/CIBERSORTxFractionsWorkflow.slurm.sh
#
# slurm argument are ignored.
#
# Job name:
#SBATCH --job-name=aedavids-pipeline.slurm.sh
#
# Partition - This is the queue it goes in:
#SBATCH --partition=long
#
# Where to send email (optional)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aedavids@ucsc.edu
#

#
# configure required resoures
#
# Number of nodes you need per job:
#SBATCH --nodes=1
#
# Memory needed for the jobs.  Try very hard to make this accurate.  DEFAULT = 4gb
# S-aedwip-BATCH --mem=4gb
#S-aedwip-OOM-BATCH --mem=32G # total memory per node
#SBATCH --mem=64G # total memory per node
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=1
#
# Processors per task:
# At least eight times the number of GPUs needed for nVidia RTX A5500
# ref https://giwiki.gi.ucsc.edu/index.php?title=Using_Docker_under_Slurm
# slurm cannot limit the resources that docker uses. 
# set to number of samples / CIBERSORTxFractionsWorkflow.wdl.input.json CIBERSORTxFractionsWorkflow.numSamplesInPartition
#SBATCH --cpus-per-task=32
#
# Number of GPUs, this can be in the format of "--gres=gpu:[1-8]", or "--gres=gpu:A5500:[1-8]" with the type included (optional)
#SBATCH --gres=gpu:1
#
# Standard output and error log
#SBATCH --output=serial_test_%j.log 
#
# Wall clock limit in hrs:min:sec:
# best500GTEx* timmed out
#S-AEDWIP-BATCH --time=10:02:30
#SBATCH --time=24:00:00
#

#
## Command(s) to run (example):
#

#
# output env info to make debugging easier
#
pwd; 
hostname; 
date

# print all the cli argumnents
printf "command line arguments \n"
for var in "$@"
do
    printf "argument :$var \n"
done

# parse the arguments
numberOfArguments=12
if [ $# -lt $numberOfArguments ];
    then
        printf "ERROR missing command line arguments. expected $numberOfArguments recevied $# \n"
        exit 1 # error
    fi

set -x
colData=$1
countData=$2
deseqResultsDir=$3
findModule=$4
estimatedScalingFactors=$5
outDir=$6
wdlImportZip=$7
wdlInputJSON=$8
wdlTools=$9
gitRoot="${10}" # weird $10 evaluates to $1 not the tenth argument
ciberSortSecurityToken="${11}"
ciberSortUser="${12}"

#
# any arguments after $12 should be passed as vargs to
# pipeline.slurm.sh
#
# https://unix.stackexchange.com/a/314041
vargs=""
i=13
while [ "$i" -le "$#" ]; do
  eval "arg=\${$i}"
  #printf '%s\n' "Arg $i: $arg"
  vargs="${vargs} $arg"
  i=$((i + 1))
done

printf "vargs: $vargs\n"


#
# 5/29/24
# pipeline.upstreamPipeline.py createSignatureMatrix() has an extra boolean flag argument
# if true use mean, else use mean to calculate expected count values for each gene in a category
# for backwards compatiblity the default value is False and we use an environmental variable
# to pass the flag
if [[ -z "${USE_MEDIAN}" ]]; then
    # True if the length of string is zero
    useMedia=""
else
  MY_SCRIPT_VARIABLE="${DEPLOY_ENV}"
  useMedian="--useMedian"
fi

printf "useMedian : ${useMedian}\n"

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
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
# run the upstream pipeline
#

printf "\n\n\nBEGIN run upstream pipeline \n"
python -m pipeline.upstreamPipeline \
            --colDataPath "${colData}" \
            --countDataPath "${countData}" \
            --deseqResultsDir "${deseqResultsDir}" \
            --estimatedScalingFactors "${estimatedScalingFactors}" \
            --findModule "${findModule}" \
            --outDir "${outDir}" \
            "${useMedian}" \
            --vargs ${vargs}

printf "END upstream pipeline  \n"

#
# run CIBERSORTx
# ref: extraCellularRNA/terra/cibersortx/wdl/CIBERSORTxFractionsWorkflow.slurm.sh
#      extraCellularRNA/terra/cibersortx/wdl/runCIBERSORTxFractionsWorkflow.sh
#
printf "\n\n\nBEGIN CIBERSORTxFractionsWorkflow.wdl\n"

signatureGenesFile=`find ${outDir} -name signatureGenes.tsv`
mixtureFile=`find ${outDir} -name mixture.txt`

# our conda environment should have the correct version of java
which java
java -version
printf "wdlTools : ${wdlTools}\n"

#java -jar "${wdlTools}/womtool-85.jar" validate CIBERSORTxFractionsWorkflow.wdl ;

wdlInputTemplate="${gitRoot}/deconvolutionAnalysis/bin/CIBERSORTxFractionsWorkflow.wdl.input.json.template"
wdlInputJSON="${outDir}/CIBERSORTxFractionsWorkflow.wdl.input.json"

#
# it is hard to replace the template tokens with file paths
# because the '/' is the separator in ed commands
# work around:
# 1. replace '/' with 'xxx'
# 2. replace token in template
# 3. replace 'xxx' with '/'
#
printf "\n\n\n configure cromwell input configuration files\n"

otd=`echo $outDir/CIBERSORTxFractionsWorkflow.wdl.output | sed -e 's/\//xxx/g'`
sg=`echo $signatureGenesFile | sed -e 's/\//xxx/g'`
mx=`echo $mixtureFile        | sed -e 's/\//xxx/g'`
cat  "${wdlInputTemplate}" | sed -e "s/TOKEN-signatureGenes.tsv/${sg}/g" \
                                -e "s/TOKEN-mixture.txt/\\${mx}/g" \
                                -e "s/TOKEN-CIBERSORTx-security/${ciberSortSecurityToken}/g" \
                                -e "s/TOKEN-CIBERSORTx-user/${ciberSortUser}/g" \
                                > t

cat t | sed  -e 's/xxx/\//g' > ${wdlInputJSON}
'rm' -rf t

printf "\n\n\n"

cat cromwellOptions.json.template | sed -e "s/TOKEN_OUTPUT_DIR/${otd}/g" \
                                    > t
                                    
cat t | sed  -e 's/xxx/\//g' > cromwellOptions.json
'rm' -f t



# make sure the docker images are on what ever phoenix node slurm assigned our job to
# phoenix will periodically remove old images
printf "\n\n\n make sure docker images are on the asigned slurm node \n"
docker pull aedavids/wdltest
docker pull aedavids/cibersortx_fractions

printf "\n\n\nBEGIN CIBERSORTxFractionsWorkflow.wdl \n"
${gitRoot}/bin/runCromwell.sh \
     -Dconfig.file="${gitRoot}/terra/wdl/cromwellDebug.conf" \
     -jar "${wdlTools}/cromwell-85.jar" \
     run \
     --imports "${wdlImportZip}" \
     --inputs "${wdlInputJSON}" \
     --options cromwellOptions.json \
     CIBERSORTxFractionsWorkflow.wdl

printf "\n\n\nEND  CIBERSORTxFractionsWorkflow.wdl\n"

#
# calculate signature matrix peformance on mixture matrix
#
printf "\n\n\nBEGIN Evaluate CIBERSORTx results as if output from k-way classifier\n"
expectedFractions=`find ${outDir} -name expectedFractions.txt`
cromwellOutDir=`grep final_workflow_outputs_dir cromwellOptions.json | cut -d : -f 2 | cut -d , -f 1`
cromwellOutDir=`echo $cromwellOutDir | tr --delete \"`
python -m analysis.metrics -e "${expectedFractions}" \
                           -f "${cromwellOutDir}/results.txt" \
                           -o "${outDir}/metrics"

#
# run extra analytics
# 

#
# createLungDeconvolutionMetricsCSV
# runs grep to find lung releated metrics
#

# p is a hack. These scripts find reference files relative to where they are installed
# todo might be clean to have best*.sh copy the bin dir locally
p="/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/bin"
script="${p}/createLungDeconvolutionMetricsCSV.sh"
scriptOutDir="${outDir}/createLungDeconvolutionMetricsCSV.sh.out"
mkdir -p "${scriptOutDir}"
$script "${outDir}/metrics/metricsRounded.csv" > "${scriptOutDir}/lungMetrics.csv"

#
# createElifeDeconvolutionMetricsCSV
# runs grep 
#
scriptOutDir="${outDir}/createElifeDeconvolutionMetricsCSV.sh.out"
mkdir -p "${scriptOutDir}"
script="${p}/createElifeDeconvolutionMetricsCSV.sh "
$script "${outDir}/metrics/metricsRounded.csv" > "${scriptOutDir}/elifeMetrics.csv"

#
# lungExploreClassificationErrors
# calculates FP, FN, and set of shared genes for tempest dataset (Lung) tissue type
#
intersectionDict="${outDir}/upsetPlot.out/*intersection.dict"
scriptOutDir="${outDir}/lungExploreClassificationErrors.sh.out"
mkdir -p "${scriptOutDir}"
script="${p}/lungExploreClassificationErrors.sh"
$script "${intersectionDict}" "${cromwellOutDir}/results.txt" "${expectedFractions}" "$scriptOutDir"


#
# flatten the confusion matrix and save as csv
#
fmeOutDir="${outDir}/findMisclassificationErrors.out"
mkdir -p "$fmeOutDir"
cmPath="${outDir}/metrics/confusionMatrix.csv"
python -m analysis.findMisclassificationErrors \
        --threshold 0 \
        --confusionMatrixPath "${cmPath}" \
        --outDir "${fmeOutDir}"

printf "\n\n\nEND Evaluate CIBERSORTx results as if output from k-way classifier\n"

printf "END ${0}"
exit 0
