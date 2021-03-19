#!/bin/bash

# download data a batch of data exRNA.org
# aedavids@ucsc.edu

export scriptName=`basename $0`

if [ $# -ne 3 ]; then
    echo "error: usage $scriptName batchDownLoadInstructionFile rootDir phoneNumber  "
    echo "example $scriptName testDownLoadInstr.aa ./kimLab/exRNA.org 7508622639@txt.att.net"
    echo "    testDownLoadInstr.aa is a batch of download parameters. see README.md for more info"
	echo "    all files and directory will be created under ./kimLab/exRNA.org"
    echo "    7508622639@txt.att.net is a phone number to send 'data is up; txt msg"
    echo ""
    exit 1
fi

# track arguments and environment setting to make more reproducable
echo cli: $0 $@

# use export else variable will not be defined by inline script embedded in nohup
export batchInstructionsFile="$1"
export rootDir="$2"
export phoneNumber=$3
export dateStamp=`dateStamp.sh`
export biFile=`basename ${batchInstructionsFile}`

export logDir="./logs"
mkdir -p ${logDir}

export scriptLog="${logDir}/$scriptName.${biFile}.${dateStamp}.out"


#set -x # turn debug trace on
#set +x # turn debug trace off

# create a session 
# this will ensure jobs continue to run even after we
# log out. It will also make easier to kill all the children
# if we need to restart for some reason

# mac does not have setsid, just use sh to test
setsid sh -c 'set -x; parseInstructionsAndDownload.sh ${batchInstructionsFile} ${rootDir} ' > $scriptLog 2>&1 &
#sh -c 'set -x; parseInstructionsAndDownload.sh ${batchInstructionsFile} ${rootDir} ' > $scriptLog 2>&1 &



