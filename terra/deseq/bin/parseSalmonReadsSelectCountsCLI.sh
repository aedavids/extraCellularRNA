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

set -x # turn debug on
# set +x # turn debug off

source ${configFile}
status=$?
if [ $status -ne 0 ]
then
    echo "ERROR unable to source ${configFile}"
    exit 1
fi

mainPython="../python/bigDataDeseq/parseSalmonReadsSelectCountsCLI.py"

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

set -x # turn debug on
# set +x #turn debug off

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

CLIENT_MODE="client"
#CLIENT_MODE="cluster"


#
# configure log levels
# https://cloud.google.com/dataproc/docs/guides/driver-output?authuser=1#configuring_logging
#        --driver-log-levels root=WARN \
echo "\n\n\n************ $0 submit job"

#
# TODO: create more batches (fewer rows) to reduce total wall clock time
# split -l 3400 sparkGTEx-trainQuantFiles.csv  sparkGTEx-trainQuantFiles.csv.
#

QUANT_FILE_CSV_ROOT="gs://${BUCKET}/GTEx/parseSalmonReadsSelectCountsCLI/quantFilesCSV/sparkGTExQuantFilesBatches"
# for i in {ab,ac,ad};  
#for i in `seq 10 18`;
#for i in `seq 19 27`;
#for i in `seq 28 36`;
#for i in `seq 37 45`;
#for i in `seq 46 54`;
for i in `seq 55 59`;
do
    LABEL_ARG_T="${LABEL_ARG},sample-set=${i}"
    echo "\n ******** ${i} ********"
    gcloud dataproc jobs submit pyspark ${mainPython} \
       --cluster=${CLUSTER_NAME} \
       --bucket ${BUCKET} \
       --region ${REGION} \
       --project ${PROJECT_ID} \
        --driver-log-levels root=WARN \
       ${LABEL_ARG_T} \
       --py-files ${extraPkg} \
       --properties=spark.submit.deployMode=${CLIENT_MODE} \
       --async \
       -- \
       --quantFilesCSV="${QUANT_FILE_CSV_ROOT}/sparkGTExQuantFiles-batch.${i}.csv" \
       --selectOnlyNumReads \
       --outputDir="gs://${BUCKET}/GTEx/spark.out"
done
      

 # --quantFilesCSV="gs://${BUCKET}/quant/spark6SampleTestQuantFiles.csv" \
 #       --outputDir="gs://${BUCKET}/quant/test6Samples"`

echo "\n\n ************"
echo "ret = ${ret}"

#(base) $ grep jobId parseSalmonReadsSelectCountsCLI.sh.tiny2-test6-select-counts.2022-01-01-t20-09.out | cut -d : -f 2
# 78fbeb1f73b8400c9b8a43f1a57c27ca

#        --quantFilesCSV="gs://${BUCKET}/GTEx/parseSalmonReadsSelectCountsCLI/quantFilesCSV/sparkGTExQuantFilesBatches/sparkQuantFile-batch.0${i}.csv" \

