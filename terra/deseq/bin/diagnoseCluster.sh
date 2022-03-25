#!/bin/sh

# https://cloud.google.com/dataproc/docs/support/diagnose-command
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 clusterName gcp-project.config"
    exit 1
fi

CLUSTER_NAME=$1
configFile=$2

source ${configFile}
status=$?
if [ $status -ne 0 ]
then
    echo "ERROR unable to source ${configFile}"
    exit 1
fi

set -x # turn debug on
# set +x # turn debug off
gcloud dataproc clusters diagnose $CLUSTER_NAME \
       --project ${PROJECT_ID} \
       --region ${REGION} 

