#!/bin/bash
# aedavids@ucsc.edu


set -x # turn trace debug log on
# set +x # turn trace debug log off

java -Dbackend.providers.Local.config.runtime-attributes='String? docker String? docker_user="$EUID"' $@

