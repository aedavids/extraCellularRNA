#!/bin/sh
# Andrew E. Davidson
# aedavids@ucsc.edu
# 10/27/2022
# ref: https://dummylabs.com/posts/2018-08-15-monit/
#


# https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca
set -euxo pipefail

# turn debug trace of
set +x

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 list of docker containers names"
    echo "use 'docker ps' to gets a list of names"
    echo "use setsid to continue running after you log out"
    echo " example: setsid sh -c '$0 eager_trucks' > $0.out 2>&1 &"
    exit 1
fi

listOfContainers=$@

# dateStamp example: 2019-12-09-23.01.43-UTC
startTime=`date "+%Y-%m-%d-%H.%M.%S-%Z%n"`


#list_of_containers="homeassistant addon_a0d7b954_appdaemon3 hassio_supervisor"
while [ 1 ]
do
    
    containers=`docker ps -f status=running --format "{{.Names}}"`
    for container in $listOfContainers
    do
        if echo $containers |grep -q $container
        then  echo "$container online " > /dev/null
        else #echo "$container offline"
             #echo "$container is offline" | mail -s "$container offline" aedavids@ucsc.edu
             endTime=`date "+%Y-%m-%d-%H.%M.%S-%Z%n"`
             body="$container is offline.\n monitoring start: $startTime \n end: $endTime \n"
             printf "$body" | mail -s "$container offline" aedavids@ucsc.edu
             exit 1
        fi

    done
    
    #sleep 30
    sleep `expr 60 \* 15`

done

exit 0
