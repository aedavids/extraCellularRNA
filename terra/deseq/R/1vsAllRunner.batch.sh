#!/bin/bash

#
# Andrew E. Davidson, aedavids@ucsc.edu
# 1/25/2022
#
# use to test long running R script running in a docker container
# docker exec --user rstudio practical_keldysh /home/rstudio/extraCellularRNA/terra/deseq/R/dockerExecTest.sh

set -x # turn debug on
whoami
cd /home/rstudio/extraCellularRNA/terra/deseq/R
export PATH=".:${PATH}"
1vsAllRunner.sh 2>&1 > 1vsAllRunner.sh.out2
exitStatus=$?
echo "exit status = $exitStatus"
exit $exitStatus

