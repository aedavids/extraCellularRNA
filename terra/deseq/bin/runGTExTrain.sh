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
createGTExTrainingCluster.sh ${CLUSTER_NAME} ./data/GTEx-training.config 2>&1 | tee createGTExTrainingCluster.sh.${CLUSTER_NAME}.out


echo "\n\n\n************ submit job"
submitSparkGTexTrain.sh ${CLUSTER_NAME} ./data/GTEx-training.config name=gtex-${CLUSTER_NAME} 2>&1 | tee submitSparkGTexTrain.sh.${CLUSTER_NAME}.out
