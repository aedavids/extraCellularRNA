#!/bin/sh

#
# ref:
# https://cloud.google.com/sdk/gcloud/reference/dataproc/jobs/submit/spark
# https://spark.apache.org/docs/latest/configuration.html#available-properties
#

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 clusterName gcp-project.config labels"
    echo "labels are optional"
    echo "example label 'env=prod,customer=acme"
    exit 1
fi

CLUSTER_NAME=$1
configFile=$2
labels=$3

set -x # turn debug on
# set +x # turn debug off

source ${configFile}
status=$?
if [ $status -ne 0 ]
then
    echo "ERROR unable to source ${configFile}"
    exit 1
fi

mainPython="../python/bigDataDeseq/preprocessData.py"

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

#CLIENT_MODE="client"
CLIENT_MODE="cluster"

# originally we did not set dirver memory. looks like default was 3840m
# we get taskMemoryManger unable to allocate 16777216 bytes when i == 2 or 3


#
# configure log levels
# https://cloud.google.com/dataproc/docs/guides/driver-output?authuser=1#configuring_logging
#        --driver-log-levels root=WARN \

# ,spark:spark.driver.memory=6g

echo "\n\n\n************ $0 submit job"
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
       --mappingCSV="gs://${BUCKET}/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv" \
       --quantFilesCSV="gs://${BUCKET}/quant/spark6SampleTestQuantFiles.csv" \
       --outputDir="gs://${BUCKET}/quant/spark6SampleTest.out"

# https://stackoverflow.com/questions/62386571/dataproc-cluster-propertycore-memory-and-memoryoverhead-setting


# https://stackoverflow.com/questions/66312936/where-do-i-configure-spark-executors-and-executor-memory-of-a-spark-job-in-a-dat
#  Feb 22 '21 at 9:58
# gcloud dataproc jobs submit pyspark path_main.py --cluster=$CLUSTER_NAME \
# --region=$REGION \
# --properties="spark.submit.deployMode"="cluster",\
# "spark.dynamicAllocation.enabled"="true",\
# "spark.shuffle.service.enabled"="true",\
# "spark.executor.memory"="15g",\
# "spark.driver.memory"="16g",\
# "spark.executor.cores"="5"


#echo "\n\n\n************ $0 submit job"
echo run gcloud dataproc jobs describe job_id
