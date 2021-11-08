#!/bin/bash
# aedavids@ucsc.edu


set -x # turn trace debug log on
# set +x # turn trace debug log off
AEDWIP

cromwell was installed by conda. I do not think this is needed anymore

java -Dbackend.providers.Local.config.runtime-attributes='String? docker String? docker_user="$EUID"' $@

