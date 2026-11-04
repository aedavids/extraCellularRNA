#!/bin/sh
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/29/2020
#

if [ $# -ne 1 ]; then
	scriptName=`basename $0`
	echo "error: usage $scriptName DESeq2 output data file"
	echo "example: $ $scriptName /public/groups/kimlab/kras.ipsc/day.7.de.seq.csv"
    echo "missing argument"
    exit 1
fi

inputFile=$1

rootDir=`git rev-parse --show-toplevel`

set -x # turn debug on
#set +x # turn debug off

# create output file n

# s#^/## If the first character (^/) is a /, remove it
# ; separates multiple sed commands
# s#/#.#g Replace all remaining / with .
# use # as the delimitor instead of / so we do not need to escape / 
path=`echo $inputFile | sed 's#^/##; s#/#.#g'`


#[[ -d img ]] || mkdir -p img
mkdir -p img

outputFile="img/${path}.png"

python -m src.bme263DataVis.volcanoPlot \
		-i ${inputFile} \
		-o ${outputFile}


