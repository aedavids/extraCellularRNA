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
# $ gsutil du -sh  gs://${BUCKET}/quant/*.quant.sf.gz
# 835.55 GiB   gs://anvil-gtex-v8-hg38-edu-ucsc-kim-lab-spark/quant/*.quant.sf.gz
#
# gzip 72% max size reduction
#  https://www.pingdom.com/blog/can-gzip-compression-really-improve-web-performance/#:~:text=Compression%20is%20a%20CPU%2Dintensive,but%20at%20a%20lower%20speed.
#
# k = compression ration = 0.72
# u = uncompressed size
# c = compressed size 835.55
#
# c = (u - c) / u
# u = c / (1 -k)
# u = 2984.10 GiB
# u = 2984.10 / 1024 ~ 3 TB
# 

#
# machine type m1-megamem-96 has 96 vcpu and 1.4 TB of memory
#        --worker-machine-type m1-megamem-96 \
#

if [ ! -z ${notebook} ]
then
    JUPYTER_NOTEBOOK="--enable-component-gateway"
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


TINY="--master-machine-type n1-standard-4 \
      --master-boot-disk-size 500 \
      --num-workers 2 \
      --worker-machine-type n1-standard-4 \
      --worker-boot-disk-size 500 \
      --num-worker-local-ssds 1 \
      "

#
# m1-megamem-96 vcpu=96 memory 1.4 tb.
# estimated cost for 2 machines, in Iowa $7.46/hr
#
MEGA_MEM="--master-machine-type n1-standard-4 \
          --master-boot-disk-size 500 \
          --num-master-local-ssds 1 \
          --num-workers 2 \
          --worker-machine-type m1-megamem-96 \
          --worker-boot-disk-size 500 \
          --num-worker-local-ssds 1 \
          "

echo "create TINY cluster"; MACHINE_TYPE="${TINY}"
# echo "create MEGA_MEM cluster"; MACHINE_TYPE="${MEGA_MEM}"

echo "\n\n\n*********** $0 create cluster"
gcloud dataproc clusters create ${CLUSTER_NAME} \
       \
       --bucket ${BUCKET} \
       --region ${REGION} \
       --zone ${ZONE} \
       ${MACHINE_TYPE} \
       --image-version 2.0-ubuntu18 \
       --properties ${LOG_PROPS} \
       --max-idle 10m \
       --project ${PROJECT_ID} \
        ${JUPYTER_NOTEBOOK}

echo "\n\n\n $0 ************"

if [ ! -z ${notebook} ]
then
    #
    # print URL to jupyter notebook server
    #
    gcloud dataproc clusters describe ${CLUSTER_NAME} --region=${REGION} --project=${PROJECT_ID} | grep jupyter
fi


