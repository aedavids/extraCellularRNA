#
# Andrew E. Davidson, aedavids@ucsc.edu
# 7/28/2020
#

# 
# starts a detached docker container
# when you are done with the container
# 1) find the container id
#   docker ps 
# 2) stop and remove the container
#   docker rm -f container_id 
#

set -x # turn debug on
# set + x # turn debug off

#PORT=`findUnusedPort.sh`
#PORT=8755
#HOST_PORT=875

HOST_PORT=`findUnusedPort.sh`
echo "ssh tunnel port number: " $HOST_PORT
CONTAINER_PORT=8787
USER_ID=`id -u`

#IMG='rocker/rstudio:4.0.0-ubuntu18.04'
#IMG='aedavids/ggplot2'
#IMG='rocker/rstudio:3.5.0'
#IMG='aedavids/biocworkshops' can not install Desq
#IMG='bioconductor/bioconductor_docker:devel'
#IMG='aedavids/biocworkshop2018desq2'
#IMG='aedavids/extra_cellular_rna'
IMG='aedavids/extra_cellular_rna_2_01' # starts rstudio-server
#IMG='aedavids/edu_ucsc_kim_lab-1vsall_1.0' # production version, support for DESeq, rstudio-server was removed

# docker arguments
# -d  --detach Run container in background and print container ID
# -rm Automatically remove the container when it exits
# --publish -p Publish a container's port(s) to the host
#	-p 127.0.0.1:80:8080/tcp
#	This binds port 8080 of the container to TCP port 80 on 127.0.0.1 of the host machine.
# --publish-all , -P Publish all exposed ports to random ports
# --read-only
# --volume Bind mount a volume
# --workdir , -w Working directory inside the container
#docker run --rm -p 127.0.0.1:${PORT}:8787 -e DISABLE_AUTH=true

set -x # turn debug on
# set +x # turn debug off



docker run --rm \
       --detach \
       --publish 127.0.0.1:${HOST_PORT}:${CONTAINER_PORT}/tcp \
       -e DISABLE_AUTH=true \
        -e USER=rstudio \
	-e USERID=${USER_ID} \
	-e PASSWORD=ggg \
	-v /private/home/${USER}:/home/rstudio \
	-v /private/groups/kimlab:/home/kimlab \
        -v /scratch/aedavids:/scratch/aedavids \
        ${IMG}

#set -x # turn debug on
set +x # turn debug off

# get the output col headers
echo "docker ps | grep ${IMG}"
docker ps | head -1

# find our images id
docker ps | grep ${IMG}
