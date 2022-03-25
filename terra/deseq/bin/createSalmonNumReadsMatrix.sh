#!/bin/sh

if [ "$#" -lt 2 ]; then
    echo "ERROR"
    echo "Usage: $0 SalmonNumReadsTransposeMatrix.csv.txt outputFile"
    echo " a file created by createSalmonNumReadsTransposeMatrix.sh"
    exit 1
fi

m=$1
outputFile=${2}

# find path to python directory
scriptPath=`which $0`
installDir=`dirname ${scriptPath}`
pythonDir="${installDir}/../python"

set -x # turn debug on
# set +x # turn debug off
#export PYTHONPATH="${PYTHONPATH}:${pythonDir}/bigDataDeseq/"
#python transpose.py $m $outputFile

# -u do not buffer stdout or stderr. makes debug easier
python -u ${pythonDir}/bigDataDeseq/transpose.py -t $m -o $outputFile

