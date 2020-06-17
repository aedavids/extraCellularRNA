#!/bin/sh
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/29/2020
#

set -x # turn debug on
#set +x # turn debug off

find /public/groups/kimlab  -name "*de.seq*" -print |grep -v "permission denied" 

