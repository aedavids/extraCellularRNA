#!/bin/sh

if [ -z "${SPARK_HOME}" ];
then
    echo "ERROR SPARK_HOME is not define "
    exit 1
fi


set -x # turn debug on
# set +x # turn debug off

# find path to python directory
scriptPath=`which $0`
installDir=`dirname ${scriptPath}`
tmp="${installDir}/../sparkHistoryLogs"
pushd $tmp
historyLogDir=`pwd`
popd

export SPARK_HISTORY_OPTS="-Dspark.history.fs.logDirectory=file:${historyLogDir}";
$SPARK_HOME/sbin/start-history-server.sh

open http://localhost:18080


echo "*********** to stop"
echo '$SPARK_HOME/sbin/stop-history-server.sh'
