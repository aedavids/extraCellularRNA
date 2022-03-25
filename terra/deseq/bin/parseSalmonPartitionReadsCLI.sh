#!/bin/sh

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 clusterName gcp-project.config labels"
    echo "labels are optional"
    echo "example label 'env=prod,customer=acme"
    exit 1
fi

CLUSTER_NAME=$1
configFile=$2
labels=$3


source ${configFile}
status=$?
if [ $status -ne 0 ]
then
    echo "ERROR unable to source ${configFile}"
    exit 1
fi

set -x # turn debug on
# set +x # turn debug off

mainPython="../python/bigDataDeseq/parseSalmonPartitionReadsCLI.py"

#
# make sure zip file is up to date
#
echo "\n\n\n************ $0 update python zip file"
extraPkg="../python/bigDataDeseq.zip"
createBigDataDeseqZip.sh
status=$?
if [ $status -ne 0 ] ;
then
    echo "ERROR createBigDataDeseqZip.sh failed"
    exit 1
fi

if [ ! -f $extraPkg ] ;
then
    echo "ERROR $extraPkg does not exist. you need to run $0 from correct location "
    exit 1
fi

if [ ! -f $mainPython ] ;
then
    echo "ERROR $mainPython does not exist. you need to run $0 from correct location "
    exit 1
fi

if [ ! -z ${labels} ]
then
    # remove whitespace
    cleanLabels=`echo ${labels} | sed 's/ //g' `
    LABEL_ARG="--labels ${cleanLabels}"
fi


#
# on mustard
# $ pwd
# /private/groups/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
# (base) [aedavids@mustard gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx]$ wc -l quant.sf 
# 5387496 quant.sf
# numRows in each part = 5387496 / 40
# numRows = 134687
#numParts=40

GTEX_TRAIN_ARGUMENTS=" --numParts=40 \
       --quantFilesCSV=gs://${BUCKET}/quant/sparkGTEx-trainQuantFiles.csv \
       --outputDir=gs://${BUCKET}/quant/sparkGTEx-train.out  "

batchNumber=0
GTEX_TRAIN_BATCH_ARGUMENTS=" --numParts=40 \
       --outputDir=gs://${BUCKET}/quant/sparkGTEx-train.out \
       --quantFilesCSV=gs://${BUCKET}/quant/createDESeqData.output/GTEx-batch-${batchNumber}-QuantFiles.csv "

SAMPLE_6_ARGUMENTS=" --numParts=3 \
       --quantFilesCSV=gs://${BUCKET}/quant/spark6SampleTestQuantFiles.csv \
       --outputDir=gs://${BUCKET}/quant/test6Samples "

#
# configure where driver runs and where driver logs go
# https://cloud.google.com/dataproc/docs/guides/driver-output
#
# spark.submit.deployMode="cluster" runs driver on a worker
# driver output is listed in YARN userlogs, which can be accessed in Logging.
#   to find / view output
#   0. go to dataproc job page and find your job id
#   1. dataproce cluster tab. select your cluster detail page
#   2. click on view logs
#   3. use pull downs to select
#       a. first column select cloud dataproc job -> us-central1 -> your job id
#   4. all logs pull down should display dataproc.job.dirver and data.proc.yarn.container
#
#   new logger query example
#     resource.type="cloud_dataproc_job"
#     resource.labels.region="us-central1"
#     resource.labels.job_id="e73c8bdfc5f44c82903c3045679d01f1"
#     logName="projects/dataprocspark-328421/logs/dataproc.job.driver"
#     severity>=WARNING
#
# spark.submit.deployMode='client' runs driver on master
# driver output is  streamed the driver for viewing
#

CLIENT_MODE="client"
#CLIENT_MODE="cluster"


#
# configure log levels
# https://cloud.google.com/dataproc/docs/guides/driver-output?authuser=1#configuring_logging
#        --driver-log-levels root=WARN \
echo "\n\n\n************ $0 submit job"

#echo "using GTEX_TRAIN_ARGUMENTS"; arguments=${GTEX_TRAIN_ARGUMENTS}
#echo "using SAMPLE_6_ARGUMENTS"; arguments=${SAMPLE_6_ARGUMENTS}
echo "using GTEX_TRAIN_BATCH_ARGUMENTS batchNumber=${batchNumber}"; arguments=${GTEX_TRAIN_BATCH_ARGUMENTS} 
gcloud dataproc jobs submit pyspark ${mainPython} \
       --cluster=${CLUSTER_NAME} \
       --bucket ${BUCKET} \
       --region ${REGION} \
       --project ${PROJECT_ID} \
        --driver-log-levels root=INFO \
       ${LABEL_ARG} \
       --py-files ${extraPkg} \
       --properties=spark.submit.deployMode=${CLIENT_MODE} \
       --async \
       -- \
       ${arguments}

