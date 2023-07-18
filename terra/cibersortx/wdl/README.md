# CIBERSORTx Fractions WDL workflow overview

```
aedavids@ucsc.edu
3/19/23
```

## Goal:
Our training mixture file is to big to upload to https://cibersortx.stanford.edu/. When we ran the CIBERSORTx fractions docker on mustard it took over 3 days to complete. Our wdl workflow splits the mixture file into parts and uses the scatter/gather pattern to run CIBERSORT on the parts in parrallel 

ref: 
- ../cibersortParallelization.md 
- ../wdlTest/README.md . scatter/gather POC
-  https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather
- ?? we download a zip file from standford that contained a lot of readme files
- you need to register with standford to use the docker. they will send you a token you pass as an argument


getting help:
- https://bioinformatics.stackexchange.com/search?q=wdl
-  https://openwdl.slack.com/ssb/redirect


**Table Of Contents**
- step 1: Create docker image 
- step 2: Create a CIBERSORTx fractions wdl task

## step 1 Create docker image
The cibersortx/fractions from cibersortx.standford.edu can not be called directly from a wdl task. It assume you have mounted your /src/data and /src/datadir to  your local file system. We do not have docker control over mounding. Also the standford docker does not expose a function we can call from our docker task. 

See dockerFile.wdlCibersort

- what is the docker image end point?
  a wdl task is just a bash script. what application does the wdl task need to call?
    ```
    # slack wdl-nightmares 3/14 glennhicky
    
    $ docker pull cibersortx/fractions
Using default tag: latest
latest: Pulling from cibersortx/fractions
Digest: sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc
Status: Image is up to date for cibersortx/fractions:latest
    ```
    
    inspect dumps a lot of information. look for "Entrypoint"
    ```
    aedavids@mustard $ docker inspect cibersortx/fractions
[
    {
        "Id": "sha256:82a17fe6bf93b6a2ae98faefedf04cde66ebc4cdc2e6539a26ac77df353c768a",
        "RepoTags": [
            "cibersortx/fractions:latest"
    ```
    
    our wdl task will need to call
    ```
    "Cmd": [
                "/bin/sh",
                "-c",
                "#(nop) ",
                "ENTRYPOINT [\"./CIBERSORTxFractions\"]"
            ],
    ```
    
    ```
     "Entrypoint": [
                "./CIBERSORTxFractions"
            ],
    ```
***
## step 2: Create a CIBERSORTx fractions wdl task

<span style="color:red">add notes about using docker biuld file in ../wdlTest partition task</span>

validate wdl and run

initialize terminal
```
conda activate extraCellularRNA
export WDL_TOOLS=/private/home/aedavids/extraCellularRNA/java/bin
```

validate
```
java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar validate cibersortxTask.wdl
```

create input json template
```
java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar inputs cibersortxTask.wdl
```


run with test data
```
runCibersortxFractionsTask.sh 2>&1 | tee runCibersortxFractionsTask.sh.out
```

you should find a results file. You should also inspect stdout
```
$ find cromwell-executions -name CIBERSORTx_aedwip_label_Results.txt
cromwell-executions/cibersortxFractionsWorkflow/7a5ee882-dec7-4ef7-95b1-8ab81df79a19/call-cibersortxFractionsTask/execution/CIBERSORTx_aedwip_label_Results.txt
```
results should look something like. cibersort uses a monte-carlo simulation. we do not set the randome seed. Likley your results will be different
```
 cat cromwell-executions/cibersortxFractionsWorkflow/7a5ee882-dec7-4ef7-95b1-8ab81df79a19/call-cibersortxFractionsTask/execution/CIBERSORTx_aedwip_label_Results.txt
Mixture	T1	T2	T3	P-value	Correlation	RMSE
S1	0.99999999999999989	0.00000000000000012	0.00000000000000000	0.00000000000000000	1.00000000000000000	0.17515043817977474
S2	0.00022852600609523	0.99931443200125780	0.00045704199264684	0.00000000000000000	0.99999994770648137	0.17471418357263008
S3	0.00000000000000000	0.00016059614163199	0.99983940385836811	0.00000000000000000	0.99999998280039615	0.33942140553181471
S4	0.50000000110937914	0.49999999889062091	0.00000000000000000	0.00000000000000000	0.97332852678457504	0.34187081755819654
S5	0.49999999967803560	0.00000000060989777	0.49999999971206660	0.00000000000000000	1.00000000000000000	0.51223165638777368
S6	0.00000000000000000	0.49999999971531511	0.50000000028468483	0.00000000000000000	0.96832966378123353	0.35133759608713955
```

```
$ find cromwell-executions -name stdout
cromwell-executions/cibersortxFractionsWorkflow/7a5ee882-dec7-4ef7-95b1-8ab81df79a19/call-cibersortxFractionsTask/execution/stdout
```


# Phoenix HPC Slurm 
when we submit a slurm job, it get assigned to 1 of the nodes in the cluster. The nodes will not have access to docker images on mustard, or phoenix. The solution is to build the image on mustard and push it to docker.io
 
ref:
    - [extraCellularRNA/terra/wdl section 3.](../../wdl/README.md)
    - [Submit a Slurm Batch Job](https://giwiki.gi.ucsc.edu/index.php/Overview_of_using_Slurm#Submit_a_Slurm_Batch_Job)
    - [https://giwiki.gi.ucsc.edu/index.php/Quick_Reference_Guide](https://giwiki.gi.ucsc.edu/index.php/Quick_Reference_Guide)
    - [https://gypsum-docs.cs.umass.edu/slurm_tutorial.html](https://gypsum-docs.cs.umass.edu/slurm_tutorial.html)

We already have a docker account

## 1. copy to staging area

```
stagingArea=/private/groups/kimlab/aedavids/slurm-jobs/cibersortx/GTEx_TCGA
mkdir -p $stagingArea

cd extraCellularRNA/terra/cibersortx/wdl
cp cibersortxFractionsTask.wdl CIBERSORTxFractionsWorkflow.wdl runCIBERSORTxFractionsWorkflow.sh $stagingArea
cp CIBERSORTxFractionsWorkflow.wdl.input.json $stagingArea
cp CIBERSORTxFractionsWorkflow.slurm.sh $stagingArea

mkdir -p ${stagingArea}/wdlTest
cp  ../wdlTest/partitionDataTask.wdl ../wdlTest/mergeTask.wdl ${stagingArea}/wdlTest

```

## 2. push img built on mustard to docker.io

docker.io periodically remove old images. To build the image on mustard. there are two images

### 2.a build docker used by ../wdlTest/partitionDataTask.wdl and ../wdlTest/mergeTask.wdl
```
cd extraCellularRNA/terra/cibersortx/wdlTest/
myTag=aedavids/wdltest
buildContexDir=`pwd`
docker build --file ./dockerFile.wdlTest --tag $myTag $buildContexDir
```

### 2.b build the docker used by cibersortxFractionsTask.wdl
```
cd extraCellularRNA/terra/cibersortx/wdl/
myTag=aedavids/cibersortx_fractions
buildContextDir=`pwd`
docker build --file ./dockerFile.wdlCibersort --tag $myTag $buildContexDir
 => exporting to image                                                                0.0s
 => => exporting layers                                                               0.0s
 => => writing image sha256:0fc13a36b7c2607bbdc11464ff4fb8b2295166ea46caa0fae2050938  0.0s
 => => naming to docker.io/aedavids/cibersortx_fractions                              0.0s
```

### 2.c push Tag=aedavids/wdltest
```
 docker login
 dockerUser=aedavids
 dockerRepo=edu_ucsc_kim_lab
 localTag=aedavids/wdltest
 docker image push $localTag
```

### 2.d push Tag=aedavids/cibersortx_fractions
```
docker login
dockerUser=aedavids
dockerRepo=edu_ucsc_kim_lab
localTag=aedavids/cibersortx_fractions
docker image push $localTag
Using default tag: latest
The push refers to repository [docker.io/aedavids/cibersortx_fractions]
2a44f8ed016f: Pushed 

...

cc4590d6a718: Mounted from cibersortx/fractions 
latest: digest: sha256:e5cad92a010bcd2d3c7ee88ff5f35d0a88e5dfea0d660a987047a39bc8089179 size: 4497
```

## 3. <span style="color:red">TODO edit CIBERSORTxFractionsWorkflow.wdl.input.json</span>

docker can only mount file paths that contain letters, numbers or "_" and "-". You can use a symbolic link to work around paths with other chars
```
{
  "CIBERSORTxFractionsWorkflow.sigmatrix": "/private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design_tilda__gender_+_category-padj-0.001-lfc:2.0-n-25/ciberSort/signatureGenes.txt",            
  "CIBERSORTxFractionsWorkflow.QN": "false",
  "CIBERSORTxFractionsWorkflow.verbose": "true",
  "CIBERSORTxFractionsWorkflow.token": "3f561ab6d4cf373d11f23d8e205b4b72",
  "CIBERSORTxFractionsWorkflow.username":  "aedavids@ucsc.edu",
  "CIBERSORTxFractionsWorkflow.perm": "100",
  "CIBERSORTxFractionsWorkflow.label": "aedwip_label",

  "CIBERSORTxFractionsWorkflow.mixture":  "/private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design_tilda__gender_+_category-padj-0.001-lfc:2.0-n-25/ciberSort/GTEx_TCGA_TrainGroupby_mixture.txt",        
  "CIBERSORTxFractionsWorkflow.numSamplesInPartition": "500",
  "CIBERSORTxFractionsWorkflow.isCSV": "false"
}
```

## 4. start slurm job

clean up results from old runs
```
rm -rf cromwell* serial*.log
```

submit job to slurm. sbatch submits the job. squeue display job information

```
sbatch CIBERSORTxFractionsWorkflow.slurm.sh; squeue; sleep 10;  tail -f serial*.log
```

### 5. monitoring
get detailed info for debugging while job is running
```
scontrol show jobid -dd <jobid>

aedavids@phoenix $ scontrol show jobid -dd 251011
JobId=251011 JobName=aedavids-wdlTest
   UserId=aedavids(30108) GroupId=prismuser(600) MCS_label=N/A
   Priority=1970 Nice=0 Account=(null) QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   DerivedExitCode=0:0
   RunTime=00:00:28 TimeLimit=06:01:00 TimeMin=N/A
   SubmitTime=2023-06-28T16:46:25 EligibleTime=2023-06-28T16:46:25
   AccrueTime=2023-06-28T16:46:25
   StartTime=2023-06-28T16:46:26 EndTime=2023-06-28T22:47:26 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2023-06-28T16:46:26 Scheduler=Main
   Partition=main AllocNode:Sid=phoenix:1486129
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=phoenix-15
   BatchHost=phoenix-15
   NumNodes=1 NumCPUs=32 NumTasks=1 CPUs/Task=32 ReqB:S:C:T=0:0:*:*
   TRES=cpu=32,mem=32G,node=1,billing=32
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   JOB_GRES=(null)
     Nodes=phoenix-15 CPU_IDs=0-31 Mem=32768 GRES=
   MinCPUsNode=32 MinMemoryNode=32G MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/private/groups/kimlab/aedavids/slurm-jobs/cibersortx/GTEx_TCGA/CIBERSORTxFractionsWorkflow.slurm.sh
   WorkDir=/private/groups/kimlab/aedavids/slurm-jobs/cibersortx/GTEx_TCGA
   StdErr=/private/groups/kimlab/aedavids/slurm-jobs/cibersortx/GTEx_TCGA/serial_test_251011.log
   StdIn=/dev/null
   StdOut=/private/groups/kimlab/aedavids/slurm-jobs/cibersortx/GTEx_TCGA/serial_test_251011.log
   Power=

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
