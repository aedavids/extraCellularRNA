#!/bin/sh
# aedavids@ucsc.edu
# 3/1/23

# script makes it easy to run batch job

set -x # turn debug trace on

../../../bin/runCromwell.sh \
    -Dconfig.file=../../wdl/cromwellDebug.conf \
        -jar ${WDL_TOOLS}/cromwell-85.jar \
        run \
        --inputs ./PAAD.vs.pancrease.1vsAllTask.input.json \
        ../../wdl/1vsAllTask.wdl

