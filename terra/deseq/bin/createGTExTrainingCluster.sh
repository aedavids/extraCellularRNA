#!/bin/sh

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 clusterName gcp-project.config"
    echo " to enable juypter notebook: Usage: $0 clusterName gcp-project.config notebook"
    exit 1
fi

CLUSTER_NAME=$1
configFile=$2
notebook=$3

set -x # turn debug on
# set +x # turn debug off

source ${configFile}
status=$?
if [ $status -ne 0 ]
then
    echo "ERROR unable to source ${configFile}"
    exit 1
fi

set -x # turn debug on
# set +x # turn debug off

#
# configure overview
# assume all data is stored in GCS bucket
# add autoscaling policy where 50% workers are preemptible
# configure for fault tollerance (need because task are expected to fail because of  preemption)
#
# num-workers = 2 & num-secondary-workers=2
#   YARN cores = 176
#   YARN memory = 528 GB
#

#
# sizing estimation
#
# ssh mustard
# cd /scratch/aedwi;
# gsutil cp gs://anvil-gtex-v8-hg38-edu-ucsc-kim-lab-spark/quant/GTEX-13RTK-1826-SM-5S2P3.quant.sf.gz .
# $ ll GTEX-13RTK-1826-SM-5S2P3.quant.sf.gz 
# -rw-r--r-- 1 aedavids prismuser 81M Jan  3 12:59 GTEX-13RTK-1826-SM-5S2P3.quant.sf.gz
#
# $ gzip -d ./GTEX-13RTK-1826-SM-5S2P3.quant.sf.gz
# ls -l ./GTEX-13RTK-1826-SM-5S2P3.quant.sf 
# -rw-r--r-- 1 aedavids prismuser 604M Jan  3 13:04 ./GTEX-13RTK-1826-SM-5S2P3.quant.sf
#
# see notebook page 403
# estimated size of names, and number of reads for 10411 training sample < 500 Gb
# (does not include spark internal represtation overhead)
#
#
# machine type m1-megamem-96 has 96 vcpu and 1.4 TB of memory
#        --worker-machine-type m1-megamem-96 \
#


if [ ! -z ${notebook} ]
then
    # always include --enable-component-gateway so we can access spark history server
    # JUPYTER_NOTEBOOK="--enable-component-gateway"
    JUPYTER_NOTEBOOK="${JUPYTER_NOTEBOOK} --properties dataproc:jupyter.notebook.gcs.dir=${BUCKET}/notebooks"
    JUPYTER_NOTEBOOK="${JUPYTER_NOTEBOOK}  --optional-components JUPYTER "
fi

# https://cloud.google.com/dataproc/docs/concepts/configuring-clusters/cluster-properties#file-prefixed_properties_table

#
# autoscaling and preemptible VM properties
# make spark more fault tolerant.
# do not use these properties while debugging or with very expensive clusters
# if your job fails (OOM) it will restart the job from the begining
# you will just thrash 9 more times
#
# properties must be comma separated with no spaces
PREMPTIBLE_PROPS="yarn:resouremanager.am.max-attempts=10"
PREMPTIBLE_PROPS="${PREMPTIBLE_PROPS},spark:spark.task.maxFailures=10"
PREMPTIBLE_PROPS="${PREMPTIBLE_PROPS},spark:spark.stage.max.ConsecutiveAttempts=10"


#
# Logging properties
# if not set I think logs will be deleted when cluster is deacivated
# see logging section in gcpSparkDebugNotes.md
#

# properties must be comma separated with no spaces
LOG_PROPS="dataproc:dataproc.logging.stackdriver.job.driver.enable=true"
LOG_PROPS=${LOG_PROPS}",dataproc:dataproc.logging.stackdriver.enable=true"
LOG_PROPS=${LOG_PROPS}",dataproc:dataproc.logging.stackdriver.job.yarn.container.enable=true"
LOG_PROPS=${LOG_PROPS}",dataproc:jobs.file-backed-output.enable=true"
# can we configure driver when we submit? I assume works are configured  correctly by default
#LOG_PROPS=${LOG_PROPS}",spark:spark-log4j=AEDWIP-filePath"

       # --autoscaling-policy autoScalePolicy-bucketStorage \

       # --master-machine-type n1-standard-4 \
       # --master-boot-disk-size 500 \
       # --num-master-local-ssds 2 \
       # --num-workers 8 \
       # --worker-machine-type n1-standard-8 \
       # --worker-boot-disk-size 500 \
       # --num-worker-local-ssds 4 \
       # --num-secondary-workers 8 \
       # --secondary-worker-boot-disk-size 500 \
       # --num-secondary-worker-local-ssds 4 \


#
# 8  YARN cores
# 24 GB YARN memory
#
TINY="--master-machine-type n1-standard-4 \
      --master-boot-disk-size 500 \
      --num-workers 2 \
      --worker-machine-type n1-standard-4 \
      --worker-boot-disk-size 500 \
      --num-worker-local-ssds 1 \
      "

#
# 32  YARN cores
# 96 GB YARN memory
#
TINY8W="--master-machine-type n1-standard-4 \
      --master-boot-disk-size 500 \
      --num-workers 8 \
      --worker-machine-type n1-standard-4 \
      --worker-boot-disk-size 500 \
      --num-worker-local-ssds 1 "

#
# 64 YARN cores YARN memory 192 GB
#
TINY16W="--master-machine-type n1-standard-4 \
      --master-boot-disk-size 500 \
      --num-workers 16 \
      --worker-machine-type n1-standard-4 \
      --worker-boot-disk-size 500 \
      --num-worker-local-ssds 1 "

#
# m1-megamem-96 vcpu=96 memory 1.4 tb.
# estimated cost for 2 machines, in Iowa $7.46/hr
# GTExTrainNumReadsMatrix.tsv is 213G on disk
#
# MEGA_MEM="--master-machine-type n1-standard-4 \
MEGA_MEM="--master-machine-type n1-standard-8 \
          --master-boot-disk-size 500 \
          --num-master-local-ssds 1 \
          --num-workers 2 \
          --worker-machine-type m1-megamem-96 \
          --worker-boot-disk-size 500 \
          --num-worker-local-ssds 2 \
          "

echo "\n\n\n*********** $0 create cluster"

#echo "create TINY 2 worker cluster 8 YARN cores YARN memory 24 GB"; MACHINE_TYPE="${TINY}"
#echo "create TINY 8 worker cluster 32 cores 96 GB memory"; MACHINE_TYPE="${TINY8W}"
#echo "create TINY 16 worker cluster 64 YARN cores YARN memory 192 GB"; MACHINE_TYPE="${TINY16W}"
echo "create MEGA_MEM 2 worker cluster ??vcpu=192 memory 2.8?? tb."; MACHINE_TYPE="${MEGA_MEM}"


# spark:spark.driver.memory=20g
# ,spark:spark.sql.autoBroadcastJoinThreshold=-1,spark:spark.sql.shuffle.partitions=576

# aedwip parmeterize properties for machine type ie 576
gcloud dataproc clusters create ${CLUSTER_NAME} \
       \
       --bucket ${BUCKET} \
       --region ${REGION} \
       --zone ${ZONE} \
       ${MACHINE_TYPE} \
       --image-version 2.0-ubuntu18 \
       --properties ${LOG_PROPS},spark:spark.sql.autoBroadcastJoinThreshold=-1,spark:spark.sql.shuffle.partitions=576,dataproc:dataproc.scheduler.max-cpu-load=1.0,dataproc:dataproc.scheduler.max-memory-used=1.0 \
       --max-idle 10m \
       --project ${PROJECT_ID} \
       --enable-component-gateway \
        ${JUPYTER_NOTEBOOK}

echo "\n\n\n $0 ************"
gcloud dataproc clusters describe ${CLUSTER_NAME} --region=${REGION} --project=${PROJECT_ID} 

if [ ! -z ${notebook} ]
then
    #
    # print URL to jupyter notebook server
    #
    echo "\n\n\n ********** jupyter URLS"
    gcloud dataproc clusters describe ${CLUSTER_NAME} --region=${REGION} --project=${PROJECT_ID} | grep jupyter
fi


