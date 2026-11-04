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


# see bash man page "SHELL BUILTIN COMMANDS" for details
# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
set -euxo pipefail
# set -e Exit immediately if a command in pipeline returns non-zero exit status
# set -u Treat unset variables and parameters as  an  error
# set -x turn debug trace on. output goes to stderr, normal output goes to stdout
# set -o print value of current options


#PORT=`findUnusedPort.sh`
#PORT=8755
#HOST_PORT=875

HOST_PORT=`findUnusedPort.sh`
echo "ssh tunnel port number: " $HOST_PORT
CONTAINER_PORT=8787
USER_ID=`id -u`   # 600 needed so we can write to our home dir
#GROUP_ID=`id -g`  # 
KIMLAB_GID=614    # we need to set group id to access /private/groups/kimlab

#IMG='rocker/rstudio:4.0.0-ubuntu18.04'
#IMG='aedavids/ggplot2'
#IMG='rocker/rstudio:3.5.0'
#IMG='aedavids/biocworkshops' can not install Desq
#IMG='bioconductor/bioconductor_docker:devel'

# dockerFile.exRNA_diease_biomarkers
IMG='aedavids/ex-rna_diease_biomarkers-v1.0'

#IMG='aedavids/biocworkshop2018desq2'

# dockerFile.extra_cellular_RNA
#IMG='aedavids/extra_cellular_rna' 

# starts rstudio-server
#IMG='aedavids/extra_cellular_rna_2_01'

#IMG='aedavids/edu_ucsc_kim_lab-1vsall_1.0' # production version, support for DESeq, rstudio-server was removed dockerFile.1vsAll

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


#
# 9/29/20
# IMG='aedavids/ex-rna_diease_biomarkers-v1.0' is very slow to start
# 'docker logs container_id' reports a lot of chmod permission problems
# you can ignore these error
# it just takes along time before you can connect to rstudio
#

# this will run rstudio but you can not access /home/kimlab/data
# docker run --rm \
#        --detach \
#        --publish 127.0.0.1:${HOST_PORT}:${CONTAINER_PORT}/tcp \
#        -e DISABLE_AUTH=true \
#         -e USER=rstudio \
# 	-e USERID=${USER_ID} \
# 	-e PASSWORD=bioc \
#         -v /private/home/${USER}:/home/rstudio \
#         -v /private/groups/kimlab:/home/kimlab \
#         -v /scratch/aedavids:/scratch/aedavids \
#         ${IMG}


# unable to connect browser to contain. I think rstudo does not start
# docker run --rm \
#        --detach \
#        --publish 127.0.0.1:${HOST_PORT}:${CONTAINER_PORT}/tcp \
#        --user ${USER_ID}:${GROUP_ID} \
#        -e DISABLE_AUTH=true \
# 	-e PASSWORD=bioc \
#         -v /private/home/${USER}:/home/rstudio \
#         -v /private/groups/kimlab:/home/kimlab \
#         -v /scratch/aedavids:/scratch/aedavids \
#         ${IMG}

# this does not work
# docker run --rm \
#        --detach \
#        --publish 127.0.0.1:${HOST_PORT}:${CONTAINER_PORT}/tcp \
#        --user ${USER_ID}:${GROUP_ID} \
#        -e DISABLE_AUTH=true \
#         -e USER=rstudio \
# 	-e USERID=${USER_ID} \
# 	-e PASSWORD=bioc \
#         -v /private/home/${USER}:/home/rstudio \
#         -v /private/groups/kimlab:/home/kimlab \
#         -v /scratch/aedavids:/scratch/aedavids \
#         ${IMG}


#
# usally you would set the docker user id and group id using
# --user ${USER_ID}:${GROUP_ID}
# this does not work with rocker/rstudio    
# you will not be able to connect to the container from your browser
# 
docker run --rm \
       --detach \
       --publish 127.0.0.1:${HOST_PORT}:${CONTAINER_PORT}/tcp \
       -e DISABLE_AUTH=true \
        -e USER=rstudio \
	-e USERID=${USER_ID} \
        -e GROUPID=${KIMLAB_GID} \
	-e PASSWORD=bioc \
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

printf "\nIt make take between 5 and 10 mins before you will be able to connect your browser\n"
