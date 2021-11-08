#!/bin/bash

trivial example of how to upload files. 
if [ "$#" -ne 1 ]; then
    echo "usage: ./upload.sh bucket-name"
    exit 1
fi

set -x

aediwp I think this get what ever the current project on gcp console is
PROJECT=$(gcloud config get-value project)
BUCKET=$1
INSTALL=gs://${BUCKET}/install.sh

#upload install
gsutil cp install.sh $INSTALL

