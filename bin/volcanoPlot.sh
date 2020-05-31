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
path=`echo $inputFile | sed 's/\/public\/groups\/kimlab\///g'`
path=`echo $path | sed 's/\\//./g'`

[[ -d img ]] || mkdir img
outputFile="img/${path}.png"

python src/bme263DataVis/volcanoPlot.py \
		-i ${inputFile} \
		-o ${outputFile} 
		