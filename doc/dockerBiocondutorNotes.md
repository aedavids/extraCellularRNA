# Re constructing our aedavids/extra_cellular_rna docker container
see: 
- /Users/andrewdavidson/googleUCSC/kimLab/docker/gettingStarted/notes.md
- /Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/dockerNotes.md
- /Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/dockerContainerList.md


## Background:

3/8/21 for unknow reason aedavids/extra_cellular_rna no longer works. When use try and use extraCellularRNA/bin/startRStudioServer.sh the container start we can connect it using 'docker exec it containerId /bin/bash'. We can not connect to it from our web browser we get an error 'R is taking longer to start than usual'

I know this work before we started the migration to mustard. We get the same behavior if we try and start the container on either plaza, courtyard, or mustard. We tried removing the .Rdata files. We are unable to find rserver log files. If connect to the container and kill the rserver it automatically restarts. we tried using /etc/init.d/r??server??? stop. we also tried using kill -9.

using 

```
[aedavids@plaza extraCellularRNA]$ docker image inspect aedavids/extra_cellular_rna

```

I noticed that there may be a couple of problems

1. user id
   - the container image has "USERID=30078", mustard has a different value for aedavids
   - we built aedavis/extra_cellular_rna from aedavids/biocworkshop2018desq2. maybe it is old?
   ```
    "Image": "aedavids/biocworkshop2018desq2",
            "Volumes": null,
            "WorkingDir": "",
            "Entrypoint": null,
            "OnBuild": null,
            "Labels": {
                "description": "Bioconductor docker image with system dependencies to install most packages.",
                "license": "Artistic-2.0",
                "maintainer": "maintainer@bioconductor.org",
                "name": "bioconductor/bioconductor_docker",
                "org.label-schema.license": "GPL-2.0",
                "org.label-schema.vcs-url": "https://github.com/rocker-org/rocker-versioned",
                "org.label-schema.vendor": "Rocker Project",
                "url": "https://github.com/Bioconductor/bioconductor_docker",
                "vendor": "Bioconductor Project",
                "version": "3.12.12"
   ```

# settting up a new container

build a docker image using latest bioconder image with required R packages installed

ref: [Docker containers for Bioconductor](https://www.bioconductor.org/help/docker/)

## Step 1) build image from docker file

```
$ extraCellularRNA/bin
$ TAG="aedavids/extra_cellular_rna_2_01"
$ docker build --file ./dockerFile.extra_cellular_RNA --tag $TAG .
```

## Step 2) rename the old docker image


```
[aedavids@plaza extraCellularRNA]$ docker image ls |head -n 1; docker image ls |grep aedavids
REPOSITORY                         TAG                                 IMAGE ID            CREATED             SIZE
aedavids/extra_cellular_rna        latest                              b8c2a26d9690        4 months ago        5.02GB
aedavids/biocworkshop2018desq2     latest                              e3d760f202cc        6 months ago        5.01GB
[aedavids@plaza extraCellularRNA]$ 
```


```
[aedavids@plaza extraCellularRNA]$ docker image tag b8c2a26d9690  aedavids/extra_cellular_rna_broken
```

```
[aedavids@plaza extraCellularRNA]$ docker image ls |head -n 1; docker image ls |grep aedavids
REPOSITORY                           TAG                                 IMAGE ID            CREATED             SIZE
aedavids/extra_cellular_rna          latest                              b8c2a26d9690        4 months ago        5.02GB
aedavids/extra_cellular_rna_broken   latest                              b8c2a26d9690        4 months ago        5.02GB
aedavids/biocworkshop2018desq2       latest                              e3d760f202cc        6 months ago        5.01GB
```

remove old tag
```
[aedavids@plaza extraCellularRNA]$ docker rmi aedavids/extra_cellular_rna
Untagged: aedavids/extra_cellular_rna:latest

```










***
# <span style="color:red">DEPRECATED</span> rest of document is obsolete
- build images by running command in a shell and then committing the images is a bad idea
- it is hard to know what was installed
- it is hard to reproduce
- it is error prone

rstudio server keeps getting jammed up the process bellow works but is slow and not easy to reproduce

## 1. download the lastest biocondutor image and test

I do not think we need to pull. when we run the container, if not local the pull happens automatically
```
[aedavids@plaza extraCellularRNA]$ docker pull bioconductor/bioconductor_docker
Using default tag: latest
latest: Pulling from bioconductor/bioconductor_docker

Digest: sha256:37e72eded9c113cffaf36a2d3740e106b2f923c1cc50cb970f0fb222cd34167e
Status: Downloaded newer image for bioconductor/bioconductor_docker:latest
```

```
[aedavids@plaza extraCellularRNA]$ docker image inspect bioconductor/bioconductor_docker

 "BIOCONDUCTOR_DOCKER_VERSION=3.12.31",
                "BIOCONDUCTOR_VERSION=3.12"
                
```

start the container using extraCellularRNA/bin/startRServer.sh. The biocondutor quick start say use image 'bioconductor/bioconductor_docker:devel' this caused a new image to be pulled down

looks like the version is the same?
```
[aedavids@plaza extraCellularRNA]$ docker image inspect bioconductor/bioconductor_docker:devel

BIOCONDUCTOR_DOCKER_VERSION=3.13.21",
"BIOCONDUCTOR_VERSION=3.13"

```

<span style="color:red"> rstudio fires up!!!!</span>


## 2. rename the old tag

```
[aedavids@plaza extraCellularRNA]$ docker image ls |head -n 1; docker image ls |grep aedavids
REPOSITORY                         TAG                                 IMAGE ID            CREATED             SIZE
aedavids/extra_cellular_rna        latest                              b8c2a26d9690        4 months ago        5.02GB
aedavids/biocworkshop2018desq2     latest                              e3d760f202cc        6 months ago        5.01GB
[aedavids@plaza extraCellularRNA]$ 
```


```
[aedavids@plaza extraCellularRNA]$ docker image tag b8c2a26d9690  aedavids/extra_cellular_rna_broken
```

```
[aedavids@plaza extraCellularRNA]$ docker image ls |head -n 1; docker image ls |grep aedavids
REPOSITORY                           TAG                                 IMAGE ID            CREATED             SIZE
aedavids/extra_cellular_rna          latest                              b8c2a26d9690        4 months ago        5.02GB
aedavids/extra_cellular_rna_broken   latest                              b8c2a26d9690        4 months ago        5.02GB
aedavids/biocworkshop2018desq2       latest                              e3d760f202cc        6 months ago        5.01GB
```

remove old tag
```
[aedavids@plaza extraCellularRNA]$ docker rmi aedavids/extra_cellular_rna
Untagged: aedavids/extra_cellular_rna:latest

```

ioconductor version 3.13 (BiocManager 1.30.10), ?BiocManager::install for help
Error in loadNamespace(x) : there is no package called ‘rmarkdown’
08 Mar 2021 20:31:55 [rsession-aedavids] ERROR r error 4 (R code execution error) [errormsg: Error in loadNamespace(x) : there is no package called ‘rmarkdown’
]; OCCURRED AT rstudio::core::Error rstudio::r::exec::{anonymous}::evaluateExpressionsUnsafe(SEXP, SEXP, SEXPREC**, rstudio::r::sexp::Protect*, rstudio::r::exec::{anonymous}::EvalType) src/cpp/r/RExec.cpp:186; LOGGED FROM: void rstudio::session::modules::rmarkdown::notebook::{anonymous}::onDocAdded(const string&) src/cpp/session/modules/rmarkdown/NotebookCache.cpp:340


## 3. install missing R packages

start the container. Use extraCellularRNA/bin/startRStudioServer.sh make sure to set IMG='bioconductor/bioconductor_docker:devel'

Install DESeq2. enter the following in the RStudio console
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
```

other R packages
```
install.packages("tidyr")
install.packages("readr")
install.packages("dplyr")
install.packages("rjson")
install.packages("rmarkdown")
install.packages("tximport")
```

## 4. save modifications

installing DESeq2 causes lot of files to compile. I do not think the image is portable
```
[aedavids@plaza ~]$ docker commit 82182783423e aedavids/extra_cellular_rna_public
sha256:f5bf579e2aa4b8b3f4deeaaa06556d4930b9a3776a5887909df1a1be1618d336
```

```
docker image ls |head -n 1; docker image ls |grep aedavids
REPOSITORY                           TAG                                 IMAGE ID            CREATED              SIZE
aedavids/extra_cellular_rna_public   latest                              f5bf579e2aa4        About a minute ago   4.38GB
aedavids/extra_cellular_rna_broken   latest                              b8c2a26d9690        4 months ago         5.02GB
aedavids/biocworkshop2018desq2       latest                              e3d760f202cc        6 months ago         5.01GB
```
