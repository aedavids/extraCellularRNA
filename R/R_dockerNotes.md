# R Docker Notes  7/29/2020
goal create a reproducable R environment

ref:
- [https://github.com/rreggiar/docker.utils/tree/plaza](https://github.com/rreggiar/docker.utils/tree/plaza)
- [https://hub.docker.com/r/rocker/rstudio](https://hub.docker.com/r/rocker/rstudio)

## Quick start
On couryard
```
$ /public/home/aedavids/extraCellularRNA/bin/startRStudioServer.sh
```

on Mac:

startRStudioServer.sh will print the port number to use to create our tunnel. assume it is 12345
```
$ sshTunnel.sh 12345 courtyard
```
now open http://localhost:12345 . You may need to wait a little while for the container to be come avalible

## original install
1) download the initial, parent, base container

```
$ docker pull rocker/rstudio:4.0.0-ubuntu18.04
```

2) use the container. See the quick start section above

3) use RStudio to install some packages

4) use commit so that installed packges become a durable part of a new container image
7f1700ff12b8 is the uuid for the running container. we give it a tag name aedavids/ggplot2
```
$ docker commit --author aedavids --pause \
    --message 'base image rocker/rstudio:4.0.0-ubuntu18.04. test install of ggplot is persitent' \
    7f1700ff12b8 
    aedavids/ggplot2
```

5) modify the startRStudioServer.sh to use the new image tag name

## getting version and history image for our images

```
(base) [aedavids@courtyard bin]$ docker history aedavids/ggplot2
IMAGE               CREATED             CREATED BY                                      SIZE                COMMENT
6d68a1f9595a        5 minutes ago       /init                                           98.7MB              base image rocker/rstudio:4.0.0-ubuntu18.04. test install of ggplot is persitent
6cea0b6cd0c2        8 days ago          /bin/sh -c #(nop)  CMD ["/init"]                0B                  
<missing>           8 days ago          /bin/sh -c #(nop)  EXPOSE 8787                  0B                  
<missing>           8 days ago          /bin/sh -c /rocker_scripts/install_pandoc.sh    577kB               
<missing>           8 days ago          /bin/sh -c /rocker_scripts/install_rstudio.sh   1.22GB              
<missing>           8 days ago          /bin/sh -c #(nop)  ENV PATH=/usr/lib/rstudio…   0B                  
<missing>           8 days ago          /bin/sh -c #(nop)  ENV RSTUDIO_VERSION=latest   0B                  
<missing>           8 days ago          /bin/sh -c #(nop)  ENV S6_VERSION=v1.21.7.0     0B                  
<missing>           8 days ago          /bin/sh -c #(nop)  LABEL org.label-schema.li…   0B                  
<missing>           8 days ago          /bin/sh -c #(nop)  CMD ["R"]                    0B                  
<missing>           8 days ago          /bin/sh -c /rocker_scripts/install_R.sh         674MB               
<missing>           8 days ago          /bin/sh -c #(nop) COPY dir:0e588ce349148b3e4…   73.6kB              
<missing>           2 weeks ago         /bin/sh -c #(nop)  ENV TZ=UTC                   0B                  
<missing>           2 weeks ago         /bin/sh -c #(nop)  ENV CRAN=https://packagem…   0B                  
<missing>           4 weeks ago         /bin/sh -c #(nop)  ENV R_HOME=/usr/local/lib…   0B                  
<missing>           4 weeks ago         /bin/sh -c #(nop)  ENV LANG=en_US.UTF-8         0B                  
<missing>           4 weeks ago         /bin/sh -c #(nop)  ENV LC_ALL=en_US.UTF-8       0B                  
<missing>           4 weeks ago         /bin/sh -c #(nop)  ENV TERM=xterm               0B                  
<missing>           4 weeks ago         /bin/sh -c #(nop)  ENV R_VERSION=4.0.0          0B                  
<missing>           4 weeks ago         /bin/sh -c #(nop)  LABEL org.label-schema.li…   0B                  
<missing>           6 weeks ago         /bin/sh -c #(nop)  CMD ["/bin/bash"]            0B                  
<missing>           6 weeks ago         /bin/sh -c mkdir -p /run/systemd && echo 'do…   7B                  
<missing>           6 weeks ago         /bin/sh -c set -xe   && echo '#!/bin/sh' > /…   745B                
<missing>           6 weeks ago         /bin/sh -c [ -z "$(apt-get indextargets)" ]     987kB               
<missing>           6 weeks ago         /bin/sh -c #(nop) ADD file:1e8d02626176dc814…   63.2MB
```


find origin of aedavids/ggplot2

 docker inspect aedavids/ggplot2 returns big json file, container is the parent/base we did commit on

```
[
    {
        "Id": "sha256:6d68a1f9595a909d708abeb2ab664aec0927b9299cd62fedb15690363609ae89",
        "RepoTags": [
            "aedavids/ggplot2:latest"
        ],
        "RepoDigests": [],
        "Parent": "sha256:6cea0b6cd0c22eb1bdd354ac08c32ee3637f52ef882ae6fe517e81124b465c3f",
        "Comment": "base image rocker/rstudio:4.0.0-ubuntu18.04. test install of ggplot is persitent",
        "Created": "2020-07-30T19:50:50.471682692Z",
        "Container": "7f1700ff12b8ef003ef99506721adcc74a5bece10e564c697044a227d6b6b9c4",

...

            "Image": "rocker/rstudio:4.0.0-ubuntu18.04",
            "Volumes": null,
            "WorkingDir": "",
            "Entrypoint": null,
            "OnBuild": null,
            "Labels": {
                "maintainer": "Carl Boettiger <cboettig@ropensci.org>",
                "org.label-schema.license": "GPL-2.0",
                "org.label-schema.vcs-url": "https://github.com/rocker-org/rocker-versioned",
                "org.label-schema.vendor": "Rocker Project"
            }

```