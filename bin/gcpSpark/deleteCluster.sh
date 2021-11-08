#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0  clusterName gcp-project.config"
    exit 1
fi

clusterName=$1
configFile=$2
# projectId=$2
# REGION=us-central1
set -x

source $configFile
gcloud dataproc clusters delete $clusterName \
       --project=${PROJECT_ID} \
       --region=${REGION}
       
