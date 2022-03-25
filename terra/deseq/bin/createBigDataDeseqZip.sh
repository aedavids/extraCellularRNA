#!/bin/sh


set -x # turn debug on
# set +x # turn debug off

# find path to python directory
scriptPath=`which $0`
installDir=`dirname ${scriptPath}`
pythonDir="${installDir}/../python/bigDataDeseq"

pushd "${pythonDir}/.."
'rm' bigDataDeseq.zip
zip bigDataDeseq.zip bigDataDeseq/* 

# list content
# unzip -vl bigDataDeseq.zip

if [ ! -f bigDataDeseq.zip ] ;
then
    echo "ERROR failed to create ${pythonDir}/../bigDataDeseq.zip "
    exit 1
fi
