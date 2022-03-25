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


export PYTHONPATH="${PYTHONPATH}:${SPARK_HOME}/python:${installDir}/bigDataDeseq"

#     --mappingCSV "${installDir}/test/data/mockTxId2GeneId.csv" \
#     --quantFilesCSV "${installDir}/test/data/mockQuantFiles.csv" \
#  --executor-memory MEM

# gs://anvil_gtex_v8_hg38_edu_ucsc_kim_lab_spark/kimlab/genomes.annotations/gencode.35

# --mappingCSV /private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv \
#     --quantFilesCSV "${installDir}/test/data/quantFiles.csv" \

${SPARK_HOME}/bin/spark-submit \
             --verbose \
             --driver-memory 20G \
             --executor-memory 10G \
             bigDataDeseq/preprocessData.py \
    --mappingCSV "${installDir}/test/data/mockTxId2GeneId.csv" \
    --quantFilesCSV "${installDir}/test/data/mockQuantFiles.csv" \
    --outputDir $outputDir
