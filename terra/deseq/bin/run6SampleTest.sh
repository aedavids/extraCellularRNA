#!/bin/sh

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 clusterName gcp-project.config "
    exit 1
fi

CLUSTER_NAME=$1
ts=`dateStamp.sh`
CLUSTER_NAME="${CLUSTER_NAME}-${ts}"
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

echo "\n\n\n************ create cluster"
create6SampleTestCluster.sh ${CLUSTER_NAME} ./data/GTEx-training.config 2>&1 | tee create6SampleTestCluster.sh.${CLUSTER_NAME}.out


echo "\n\n\n************ submit job"
submitSpark6SampleTest.sh ${CLUSTER_NAME} ./data/GTEx-training.config name=sample6-${CLUSTER_NAME} 2>&1 | tee submitSpark6SampleTest.sh.${CLUSTER_NAME}.out
