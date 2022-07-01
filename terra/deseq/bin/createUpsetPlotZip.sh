#!/bin/sh

if [ "$#" -ne 2 ]; then
    echo "creates zip file with all the python files need to create upset plots"
    echo "ERROR: Usage: $0 extraCellularRNARoot zipFileName"
    exit 1
fi

set -x # turn debug trace on 
root=$1
outFile=$2

rm -rf tmp
mkdir tmp

cp -r ${root}/../unmappedReadsAnalysis/python tmp
cp -r ${root}/terra/deseq/python/plots tmp/python


cd tmp
zip -r ../${outFile} python
