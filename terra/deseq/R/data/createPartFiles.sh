#!/bin/sh


if [ -z "${SPARK_HOME}" ];
then
    echo "ERROR SPARK_HOME is not define "
    exit 1
fi

set -x

outputDir="output"
mkdir -p $outputDir

t=`which $0`
installDir=`dirname t`

pythonDir="${installDir}/../../python"
pythonDir="/private/home/aedavids/extraCellularRNA/terra/deseq/python"

export PYTHONPATH="${PYTHONPATH}:${SPARK_HOME}/python"
export PYTHONPATH="${PYTHONPATH}:${pythonDir}"
#export PYTHONPATH="${PYTHONPATH}:~aedavids/extraCellularRNA/terra/deseq/python/bigDataDeseq"

echo $PYTHONPATH

"${SPARK_HOME}/bin/spark-submit" \
    --driver-memory 10g \
    --executor-memory 10g  \
    "${pythonDir}/bigDataDeseq/bigDataDeseq.py" \
    --mappingCSV "${installDir}/gencode.v35.tx.to.gene.csv" \
    --quantFilesCSV "${installDir}/quantFiles.csv" \
    --outputDir $outputDir
