#!/bin/bash
# Andrew Davidson
# aedavids@ucsc.edu
# 12/20/23
#

numArgs=4
if [ $# -ne $numArgs ];
    then
        printf "ERROR  \n"
        printf "\texpected $numArgs arguments. received $# \n"
        printf "\targs: $@ \n"
        printf "\n"
        printf "\tusage: $0 metric booleanOperation value output from python -m analysis.metrics\n"
        printf "\tmetric can be one of the following precision,recall,f1-score,support,specificity,sensitivity,tp,fn,fp,tn\n"
        printf "\tbooleanOperation can be one of the following '==', '<=', '>=' \n"
        printf "\tvalue is a number \n"
        printf "\texample: sensitivity '>=' 0.7 best25GTEx_TCGA.sh.out/metrics/metricsRounded.csv \n"
        exit 1 # error
fi

# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euo pipefail
# set -e Exit immediately if a pipeline see shell builtin command it is more complicated
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options


metric=$1
op=$2
value=$3
metricsFile=$4


# scriptPath=`which $0`
# installDir=`dirname $scriptPath`
# grep -i `cat "${installDir}/../data/elifeGrepRegx.txt"` "${metricsFile}"
#awk -F, '{if($3>=0.7)print$1, $2, $3}' < $metricsFile

#
# map metric to column number
# becareful spacing matters in bach
#
col=999
if [ $metric == "precision" ]; then
    col=2
elif [ $metric == "recall" ]; then
    col=3
elif [ $metric == "f1-score" ]; then
    col=4
elif [ $metric == "support" ]; then
    col=5
elif [ $metric == "specificity" ]; then
    col=6
elif [ $metric == "sensitivity" ];then
    col=7
elif [ $metric == "tp" ]; then
    col=8
elif [ $metric == "fn" ]; then
    col=9
elif [ $metric == "fp" ]; then
    col=10
elif [ $metric == "tn" ]; then
    col=11
elif [ $metric == ]; then
    col=12
else 
    echo "ERROR $metric is not a valid metric"
    exit 1
fi

#echo "metric: $metric, col : $col"
# awkCmd=`printf "{if(\\$${col}${op}${value})print\\$1, \\$2}" `
awkCmd=`printf "{if(\\$${col}${op}${value})print\\$1, \\$2}" `
#echo "$awkCmd"

head -n 1 "${metricsFile}" | cut -d , -f 1,${col}
#set -x
awk -F,  "${awkCmd}" < $metricsFile