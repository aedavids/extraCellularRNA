
see /Users/andrewdavidson/googleUCSC/kimLab/docker/gettingStarted/notes.md
see /Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/dockerContainerList.md

- [https://docker-curriculum.com/ good quick start](https://docker-curriculum.com/)
- [https://www.docker.com/101-tutorial](https://www.docker.com/101-tutorial)
- [more complete tutorial](https://takacsmark.com/dockerfile-tutorial-by-example-dockerfile-best-practices-2018/)

todo: pick up at [http://localhost/tutorial/multi-container-apps/](http://localhost/tutorial/multi-container-apps/)

## start tutorial

arguments:
```
-d Run container in background and print container ID
-p Publish a container's port(s) to the host
```

```
$ docker run -d -p 80:80 docker/getting-started
$ open localhost:80
```


## building the container image
- [http://localhost/tutorial/our-application/](http://localhost/tutorial/our-application/)
1. created Dockerfile
  * FROM node:12-alpine
    + Apline is small secure linux distribution
    
2. build

   ```
   $ docker build -t getting-started .
   ```
3. start contain
   ```
   $ docker run -dp 3000:3000 getting-started
   $ open http://localhost:30000
   ```
   
   
## update app
1. change one of the source files
2. build
   ```
   $ docker build -t getting-started .
   ```

## launch/run the updated container image
1. get container id
   ```
   $ docker ps
   CONTAINER ID        IMAGE                    COMMAND                  CREATED             STATUS              PORTS                    NAMES
b13248cc1563        c4fd1939735c             "docker-entrypoint.s…"   14 minutes ago      Up 14 minutes       0.0.0.0:3000->3000/tcp   suspicious_poincare
0947502ca3d1        docker/getting-started   "/docker-entrypoint.…"   27 minutes ago      Up 27 minutes       0.0.0.0:80->80/tcp       jolly_panini
   ```

2. stop  and remove
   ```
   docker stop b13248cc1563
   docker rm b13248cc1563
   ```
   
   better way to do this
   ```
   docker rm -f container_id # will stop and remove
   ```
3. start update container image
   ```
   docker run -dp 3000:3000 getting-started
   ```
   
   
## Sharing
1. create repo on [Docker hub](https://hub.docker.com/)
2. tag our image so we can push it
   a. login docker login -u aedavids
   b. docker tag getting-started aedavids/getting-started
   c. docker push aedavids/getting-started

3. running img on new instance


## persisting our DB
1. start and ubuntu container
   - bash command, pick a random num and write to file
   - docker run -d ubuntu bash -c "shuf -i 1-10000 -n 1 -o /data.txt && tail -f /dev/null"
   - using the docker dashboard we can start an interactive shell and cat teh data.txt file

We could connect from a terminal instead of using docker dashboard
```
docker exec <container-id> cat /data.txt
```


2. create named volume
this works well if we want to store data but do not care where it is. If we care see bind mounts section bellow
```
docker volume create todo-db
```

3. run container using -v to mount volume
docker run -dp 3000:3000 -v todo-db:/etc/todos getting-started

4. more about volumns
the Mountpoint is relative to the VM not our mac.
```
docker volume inspect todo-db
[
    {
        "CreatedAt": "2020-07-27T00:21:38Z",
        "Driver": "local",
        "Labels": {},
        "Mountpoint": "/var/lib/docker/volumes/todo-db/_data",
        "Name": "todo-db",
        "Options": {},
        "Scope": "local"
    }
]
$ 
```

to find out more info on our mac

- aguments
  * -rm // Automatically remove the container when it exits
  * -i // interactive Keep STDIN open even if not attached
  * -t // Allocate a pseudo-TTY
  * -v Bind mount a volume

```
$ docker run --rm -it -v /:/vm-root alpine:edge ls -l /vm-root

```

out example
```
docker run --rm -it -v /:/vm-root alpine:edge ls -l /vm-root/var/lib/docker/volumes/
total 28
-rw-------    1 root     root         32768 Jul 27 00:19 metadata.db
drwxr-xr-x    3 root     root          4096 Jul 27 00:19 todo-db
```


## using bind mounts
allows us to control exact mount point on host.

## starting a dev mode container

arguments:
-d Run container in background and print container ID
-p Publish a container's port(s) to the host
```
docker run -dp 3000:3000 \
    -w /app -v "$(pwd):/app" \
    node:12-alpine \
    sh -c "yarn install && yarn run dev"
```

Now if we make a change to a file in  apps/ the container will automatically restart. 

when we happy with the new version of our container use build

```
docker build -t getting-started .
```

??? do we need to push our image back to the docker repo?


## ssh into a running container
[https://phase2.github.io/devtools/common-tasks/ssh-into-a-container/](https://phase2.github.io/devtools/common-tasks/ssh-into-a-container/)

```
$ docker exec -it <container name> /bin/bash 
```

## saving changes made to a running container
E.G. you install a package in R studio
[https://docs.docker.com/engine/reference/commandline/commit/](https://docs.docker.com/engine/reference/commandline/commit/)
```
$ docker commit [CONTAINER_ID] [new_image_name]
```

## saving and restoring docker images (moving images from public to private machines)
original goal, docker image was created on courtyard wanted to run it on plaza
 
 ref:  BME Notebook # 2, 10/21/20 p 11 'Saving and Restoring Docker images'
 
 for following example the source machine will be plaza. The destination/target machine will be mustard
 
 1. find the image you want to save on the source machine. 
 By convention we tag our images with our email
    ```
 [aedavids@plaza ~]$ docker images |grep aedavids
 aedavids/extra_cellular_rna      latest          b8c2a26d9690        3 months ago        5.02GB
 aedavids/biocworkshop2018desq2   latest          e3d760f202cc        6 months ago        5.01GB
    ```

2. create a director on local machine to save image
   ```
   $ mkdir -p /scratch/aedavids
   ```
   
3. create a compressed tar file back up of image.
This is slow
   ```
   [plaza ~]$ docker save b8c2a26d9690 | gzip > /scratch/aedavids/extra_cellular_rna.tar.gz
   $ ls -lh /scratch/aedavids/extra_cellular_rna.tar.gz
   -rw-r--r-- 1 aedavids giuser 2.0G Feb 12 12:32 /scratch/aedavids/extra_cellular_rna.tar.gz
   ```

4. copy image to a tempory place on the destination machine
   ```
   [aedavids@mustard tmp]$ rsync -avhz aedavids@plaza:/scratch/aedavids/extra_cellular_rna.tar.gz .
   aedavids@plaza's password: 
   receiving incremental file list
   extra_cellular_rna.tar.gz

   sent 43 bytes  received 2.03G bytes  9.50M bytes/sec
   total size is 2.05G  speedup is 1.01
   [aedavids@mustard tmp]$ ls -lh ~/tmp
   total 3.9G
   -rw-r--r-- 1 aedavids prismuser 2.0G Feb 12 12:32 extra_cellular_rna.tar.gz
 
   ```

5. change the file permision
   ```
   $ chmod 755 extra_cellular_rna.tar.gz
   ```
 
 6. restore the image on the destination machine
 you should see all the layers of your image load
    ```
    $ docker load --input extra_cellular_rna.tar.gz
    b187ff70b2e4: Loading layer [==============================================>]  65.58MB/65.58MB
    
    ...
    
    Loaded image ID: sha256:b8c2a26d9690b742cbe2a8f3008cf6f1a61f560c7e826a44250916d2e60ab790
    ```
    
7. check
 ```
 $ docker images | grep b8c2a26d9690

 <none>   <none>    b8c2a26d9690        3 months ago        5.02GB
 ```
 
 8. create a tag to make it easier to work with your image
 ```
 $ docker tag b8c2a26d9690 aedavids/extra_cellular_rna
$ docker images |grep aedavids
aedavids/extra_cellular_rna   latest  b8c2a26d9690        3 months ago        5.02GB

 ```
 
 
 9. remove all the tmp files
```
plaza $ rm -i /scratch/aedavids/*
mustard $ rm -i ~/tmp/*
```


10. get runtime stats from with the container
for unknow reason app.terra.bio does not seem to respect our runtime configuration.
We need to give the container more memory to prevent salmon from crashing. [https://www.datadoghq.com/blog/how-to-collect-docker-metrics/
](https://www.datadoghq.com/blog/how-to-collect-docker-metrics/
)

as a test I used extraCellularRNA/bin/startRStudioServer.sh to start a docker container on mustard.

ssh onto mustard then connect to the contain
```
$ docker exec -it <container name> /bin/bash
```

get memory releated stats
```
root@be45d0a184be:/# cat /sys/fs/cgroup/memory/memory.kmem.limit_in_bytes 
9223372036854771712

root@be45d0a184be:/# cat /sys/fs/cgroup/memory/memory.stat 
cache 139583488
rss 4206592
rss_huge 0
mapped_file 6574080
swap 0
pgpgin 368752
pgpgout 333647
pgfault 1484463
pgmajfault 402
inactive_anon 2228224
active_anon 1929216
inactive_file 41877504
active_file 97697792
unevictable 0
hierarchical_memory_limit 9223372036854771712
hierarchical_memsw_limit 9223372036854771712
total_cache 139583488
total_rss 4206592
total_rss_huge 0
total_mapped_file 6574080
total_swap 0
total_pgpgin 368752
total_pgpgout 333647
total_pgfault 1484463
total_pgmajfault 402
total_inactive_anon 2228224
total_active_anon 1929216
total_inactive_file 41877504
total_active_file 97697792
total_unevictable 0
root@be45d0a184be:/# 
```
