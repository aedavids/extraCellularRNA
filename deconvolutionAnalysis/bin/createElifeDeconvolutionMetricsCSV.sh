#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 11/29/23
#

if [ $# -ne 1 ];
    then
        printf "ERROR  \n"
        printf "usage: $0 output from python -m analysis.metrics\n"
        printf "example: best25GTEx_TCGA.sh.out/metrics/metricsRounded.csv \n"
        exit 1 # error
fi

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
# set -euxo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options

metricsFile=$1

head -n 1 "${metricsFile}"

scriptPath=`which $0`
installDir=`dirname $scriptPath`
grep -i `cat "${installDir}/../data/elifeGrepRegx.txt"` "${metricsFile}"
