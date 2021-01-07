#!/bin/bash

fullPathToSelf=$0
if [ $# -lt 1 ]; then
    scriptName=`basename $0`
    echo "sends txt message to phone number"
    echo "error: usage $scriptName phoneNumber ..."
    echo "missing phone number. example phoneNumber 6508622639@txt.att.net"
    exit 1
fi

phoneNumber=$1
set -x

# dateStamp example: 2019-12-09-23.01.43-UTC
dateStamp=`date "+%Y-%m-%d-%H.%M.%S-%Z%n"`
host=`hostname`
echo "$host data is up $dateStamp ${@:2} " | mail ${phoneNumber}
