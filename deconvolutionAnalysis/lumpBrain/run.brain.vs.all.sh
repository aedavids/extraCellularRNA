#!/bin/bash
# aedavids@ucsc.edu
# 5/4/24

# script makes it easy to run batch job

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
#set -euxo pipefail

# do not set -u conda.sh will fail
set -exo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

printf "\n\n\nset up Java environment \n"

#
# the conda enviroment is configured to use 
# java --version
# openjdk 11.0.1 2018-10-16 LTS
# OpenJDK Runtime Environment Zulu11.2+3 (build 11.0.1+13-LTS)
# OpenJDK 64-Bit Server VM Zulu11.2+3 (build 11.0.1+13-LTS, mixed mode)
#

#we should not need to set up JAVA_HOME
# conda.sh fails if not set because of set -u above
#export JAVA_HOME=/private/home/aedavids/miniconda3/envs/extraCellularRNA

# start conda env
condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`
source ${condaBase}/etc/profile.d/conda.sh
# set -x
conda activate extraCellularRNA

gitRoot=`git rev-parse --show-toplevel`
wdlTools="${gitRoot}/java/bin"
wdl="${gitRoot}/terra/wdl/1vsAllTask.wdl"

printf "BEGIN $0"
${gitRoot}/bin/dateStamp.sh

${gitRoot}/bin/runCromwell.sh \
     -Dconfig.file="${gitRoot}/terra/wdl/cromwellDebug.conf" \
     -jar "${wdlTools}/cromwell-85.jar" \
     run \
     --inputs brain.vs.all.input.json \
     --options brain.cromwellOption.json \
     $wdl


exitStatus=$?
${gitRoot}/bin/dateStamp.sh
printf "END $0 exitStatus : ${exitStatus}\n"
