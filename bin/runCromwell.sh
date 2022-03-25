#!/bin/bash
# aedavids@ucsc.edu

if [ $# -lt 2 ]; then
	scriptName=`basename $0`
	echo "error: usage $scriptName -Dconfig.file=cromwellConfigFile  -jar cromwell-74.jar run --inputs myTask.wdl.inputs.json myTask.wdl"
    exit 1
fi

set -x # turn trace debug log on
# set +x # turn trace debug log off

java \
    -Dbackend.providers.Local.config.runtime-attributes='String? docker String? docker_user="$EUID"' \
    $@

