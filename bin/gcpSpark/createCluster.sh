#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: ./createCluster.sh clusterName gcp-project.config"
    exit 1
fi

# https://cloud.google.com/dataproc/docs/concepts/versioning/dataproc-release-2.0

CLUSTER_NAME=$1
configFile=$2

set -x

source ${configFile}

#PROJECT= #$(gcloud config get-value project)
# BUCKET=$1
# PROJECT_ID=$2 # example dataprocspark-328421
# CLUSTER_NAME=$3 # example terra-access-test1
# REGION=us-central1
# ZONE=us-central1-c

# we use INSTALL to cause extra packages to be install on all the notes
INSTALL=gs://${BUCKET}/install.sh

#upload install
gsutil cp install.sh $INSTALL

# use pricing calculate to estimate cost
# gcloud dataproc clusters create ${CLUSTER_NAME} \
#   --properties=dataproc:dataproc.personal-auth.user=aedavids@ucsc.edu \       
#   --num-workers=2 \
#   --worker-machine-type=n1-standard \
#   --master-machine-type=n1-standard \
#   --enable-component-gateway \
#   --image-version=2.0.22-ubuntu18 \
#   --optional-components=JUPYTER \
#   --project=$PROJECT_ID \
#   --region=$REGION \
#   --initialization-actions=${INSTALL}

# gcloud config set project project-id
# gcloud auth login
#  2.0-ubuntu18
# gcloud dataproc clusters create $CLUSTER_NAME \
#        --properties=dataproc:dataproc.personal-auth.user=aedavids@ucsc.edu \
#        --enable-component-gateway \
#        --region ${REGION} \
#        --zone ${ZONE} \
#        --master-machine-type n1-standard-4 \
#        --master-boot-disk-size 500 \
#        --num-workers 2\
#        --worker-machine-type n1-standard-4 \
#        --worker-boot-disk-size 500 \
#        --image-version 2.0.22-ubuntu18 \
#        --optional-components JUPYTER \
#        --project $PROJECT_ID \
#        --initialization-actions=${INSTALL}


# https://cloud.google.com/iam/docs/service-accounts-actas

# + gcloud dataproc clusters create cluster-209f --properties=dataproc:dataproc.personal-auth.user=aedavids@ucsc.edu --enable-component-gateway --region=us-central1 --optional-components=JUPYTER --project=dataprocspark-32842
# ERROR: (gcloud.dataproc.clusters.create) PERMISSION_DENIED: Permission denied on resource project dataprocspark-32842.
# - '@type': type.googleapis.com/google.rpc.Help
#   links:
#   - description: Google developer console API key
#     url: https://console.developers.google.com/project/dataprocspark-32842/apiui/credential
# - '@type': type.googleapis.com/google.rpc.ErrorInfo
#   domain: googleapis.com
#   metadata:
#     consumer: projects/dataprocspark-32842
#     service: dataproc.googleapis.com
#   reason: CONSUMER_INVALID
# (base) $ 


#        --properties=dataproc:dataproc.personal-auth.user=aedavids@ucsc.edu \
gcloud dataproc clusters create $CLUSTER_NAME \
       --enable-component-gateway \
       --region=${REGION} \
       --optional-components=JUPYTER \
       --project=$PROJECT_ID 

