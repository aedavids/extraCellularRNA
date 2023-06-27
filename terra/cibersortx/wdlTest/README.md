# CIBERSORTx scatter/gather WDL workflow POC

```
Andrew E. Davidson
aedavids@ucsc.edu
3/21/23
```

## Overview
demonstrates how we can implement a version of cibersort that splits the mixture file into parts that can be processed in parrallel. This should reduce our runtime from days to hours. testScatterGather.wdl is the complete POC workflow. You can use runTestScatterGather.sh to run the POC workflow. Each individual task is implemented in its own wdl file. You can find a run*.sh for each task. Use the run*.sh to test individual task. 

**Table of Contents**

- step 1: explore Alpine, a simple, small docker image
- step 2: create an docker image for our POC that includes some python modules we want to call from our wdl tasks
- step 3: how to create, wdl input.json, validate, and run our wdl task
- step 4, 5, 6 repeat of step three for our createDataTask, partitionDataTask, AggregateDataTask, and mergeDataTask.
  * partiitionDataTask scatters. ie implements distributed concurrency
  * mergeDataTask gathers: waits for all the distributed jobs to complete and concatenates the results into a single data file
  

# 1. create a docker img with python scripts.
we could write everything in bash. it would be kind of trick.

Alpine is a small linux based image. This one has pandas installed
```
docker pull nickgryg/alpine-pandas
```

figure out what version of python is installed, existing users, ...

```
## start docker and create an interative shell and teminal
docker run -it nickgryg/alpine-pandas ./bin/sh
```

```
/ # which python
/usr/local/bin/python
/ # python --version
Python 3.10.4
```

looks like the only user is root
```
/ # ls /home
/ #
```

use ^d to kill sh and cause docker to terminate.

```
docker ps -a # lis all dockers
docker rm -f # kill and remove docker
```

## 2. create a python script and add it to our docker image

```
myTAG=wdltest
docker build --file ./dockerFile.wdlTest --tag $myTAG .
Sending build context to Docker daemon  31.23kB
Step 1/4 : FROM nickgryg/alpine-pandas
 ---> 7368e1c031e0
Step 2/4 : MAINTAINER aedavids@ucsc.edu
 ---> Using cache
 ---> 2cc70cdc9506
Step 3/4 : COPY ./createTestData.py /bin/
 ---> Using cache
 ---> 58150f3842f0
Step 4/4 : RUN chmod a+x /bin/createTestData.py
 ---> Running in 6aae32c3752d
Removing intermediate container 6aae32c3752d
 ---> ba676e97df9c
Successfully built ba676e97df9c
Successfully tagged aedavids/wdltest:latest
```

## 3. create a createTestDataTask.wdl file and test our docker image

set up terminal env
```
export WDL_TOOLS=/private/home/aedavids/extraCellularRNA/java/bin

# use java version in conda env extraCellularRNA
conda activate extraCellularRNA
```

validate wdl
```
java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar validate createTestDataTask.wdl
```

use womtools to generate input.json template
```
java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar inputs ./createTestDataTask.wdl > createTestDataTask.wdl.input.json
```

edit. this is what json should look like
```
cat createTestDataTask.wdl.input.json
{
  "createTestWorkflow.createTestFile.n": "100"
}
```

** list of test run scripts**
```
ls *.sh
runAggregateTask.sh*   runMergeTest.sh*         runTestScatterGather.sh*  taskSplit.sh*
runCreateTestData.sh*  runParitionDataTask.sh*  task1.sh*
```

** run our test **
do not worry about exceptions when cromwell shuts down
```
/private/home/aedavids/extraCellularRNA/bin/runCromwell.sh \
     -Dconfig.file=/private/home/aedavids/extraCellularRNA/terra/wdl/cromwellDebug.conf \
     -jar ${WDL_TOOLS}/cromwell-85.jar \
     run \
     --inputs createTestDataTask.wdl.input.json \
     createTestDataTask.wdl
```

find the results
```
 $ !find
find . -name '*.test.tsv'
./cromwell-executions/createTestWorkflow/0173a0f5-74c4-4dd6-b715-d66c96d591e7/call-createTestFile/execution/100.test.tsv
```



## 4. create trival scatter gather example
Make sure we understand basic wdl syntax

see scatterGather.wdl

```
/private/home/aedavids/extraCellularRNA/bin/runCromwell.sh \
     -Dconfig.file=/private/home/aedavids/extraCellularRNA/terra/wdl/cromwellDebug.conf \
     -jar ${WDL_TOOLS}/cromwell-85.jar \
     run \
     scatterGather.wdl
```


## 5. addition test task and scatter (parrallel distibuted processing)
Scatter works. we can import local wdl files using zip. [https://cromwell.readthedocs.io/en/stable/Imports/](https://cromwell.readthedocs.io/en/stable/Imports/)

```
 zip import.zip aggregateTask.wdl  createTestDataTask.wdl  partitionDataTask.wdl mergeTask.wdl
  adding: aggregateTask.wdl (deflated 55%)
  adding: createTestDataTask.wdl (deflated 54%)
  adding: partitionDataTask.wdl (deflated 56%)
  adding: testScatterGather.wdl (deflated 58%)
  adding: mergeTask.wdl (defalted 56%)
```

```
java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar validate testScatterGather.wdl
```


Allways create a zip file, else any changes you make to the wdl files will not take effect

<span style="color:red">use runTestScatterGather.sh</span>
```
rm -rf cromwell-* import.zip; \
zip import.zip aggregateTask.wdl  createTestDataTask.wdl  partitionDataTask.wdl mergeTask.wdl; \
/private/home/aedavids/extraCellularRNA/bin/runCromwell.sh \
     -Dconfig.file=/private/home/aedavids/extraCellularRNA/terra/wdl/cromwellDebug.conf \
     -jar ${WDL_TOOLS}/cromwell-85.jar \
     run \
     --inputs testScatterGather.wdl.input.json \
     --imports /private/home/aedavids/extraCellularRNA/terra/cibersortx/wdlTest/import.zip \
      testScatterGather.wdl
```

## 5. Phoenix HPC Slurm test
when we submit a slurm job, it get assigned to 1 of the nodes in the cluster. The nodes will not have access to docker images on mustard, or phoenix. The solution is to build the image on mustard and push it to docker.io
 
ref:
    - [extraCellularRNA/terra/wdl section 3.](../../wdl/README.md)
    - [Submit a Slurm Batch Job](https://giwiki.gi.ucsc.edu/index.php/Overview_of_using_Slurm#Submit_a_Slurm_Batch_Job)
    - [https://giwiki.gi.ucsc.edu/index.php/Quick_Reference_Guide](https://giwiki.gi.ucsc.edu/index.php/Quick_Reference_Guide)
    - [https://gypsum-docs.cs.umass.edu/slurm_tutorial.html](https://gypsum-docs.cs.umass.edu/slurm_tutorial.html)

We already have a docker account

### 5a. push img built on mustard to docker.io
```
$ docker login
$ dockerUser=aedavids
$ dockerRepo=edu_ucsc_kim_lab
$ localTag=aedavids/wdltest
$ docker image push $localTag
Using default tag: latest
The push refers to repository [docker.io/aedavids/wdltest]
62d3031b59c9: Pushed 
ce07acafa4ff: Pushed 
402cc21122d3: Pushed 
dedfcc20031d: Pushed 
2bf409592244: Pushed 
7437c2bf7b24: Pushed 
9c86362ee8c7: Pushed 
851ae0d9ae05: Pushed 
971a4a1833dd: Pushed 
b18795da5e7f: Mounted from nickgryg/alpine-pandas 
dafede43d7ef: Mounted from nickgryg/alpine-pandas 
1468e3f86def: Mounted from nickgryg/alpine-pandas 
058c2739fa73: Mounted from nickgryg/alpine-pandas 
5a4db5caa17b: Mounted from nickgryg/alpine-pandas 
8cc7586db7bd: Mounted from nickgryg/alpine-pandas 
4fc242d58285: Mounted from nickgryg/alpine-pandas 
latest: digest: sha256:179d042cc08e8a4c540ebdac0c9627e14b55e56a817559c65f0011f69fc346ab size: 3659
```

### 5b copy file to slurm frienlyd staging area
You can still access file on /private/home/aedavids how ever best practices for large files to access them from /private/groups/kimlab

```
ssh from mac to phoenix
cd /private/groups/kimlab/aedavids/slurm-jobs/cibersortx/wdlTest/

# copy files to a phoenix node friendly location

d=/private/groups/kimlab/aedavids/slurm-jobs/cibersortx/wdlTest
cp $d/*.wdl .
cp $d/dockerFile.wdlTest  .
mkdir -p wdlTools
cp /private/home/aedavids/extraCellularRNA/java/bin/*.jar wdlTools/
```

clean up results from old runs
```
rm -rf cromwell* serial*.log
```

submit job to slurm. sbatch submits the job. squeue display job information

```
sbatch wdlTest.slurm.sh; squeue; tail -f serial*.log
```

### 5c monitoring
get detailed info for debugging while job is running
```
scontrol show jobid -dd <jobid>
```

after job completes
```
sacct -j <jobid> --format=JobID,JobName,MaxRSS,Elapsed
```

example
```
sacct -j 248836 --format=JobID,ReqNodes,allocNodes,NNodes,ReqMem,MaxRSS,ReqCPUS,AllocCPUS,MinCPU,NCPUS
JobID        ReqNodes AllocNodes   NNodes     ReqMem     MaxRSS  ReqCPUS  AllocCPUS     MinCPU      NCPUS 
------------ -------- ---------- -------- ---------- ---------- -------- ---------- ---------- ---------- 
248836              1          1        1        32G                  32         32                    32 
248836.batch        1          1        1              1026812K       32         32   00:00:43         32 

```
