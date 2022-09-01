# testing wdl
```
Andrew Davidson
aedavids@ucsc.edu
```

## ref:
- [5 min intro to cromwell](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/)
- [WDL github repo](https://github.com/openwdl/wdl)
- [getting started](https://support.terra.bio/hc/en-us/sections/360007274612-WDL-Documentation)



cromwell and womtool are used to write WDL scripts for terra. [download](https://github.com/broadinstitute/cromwell/releases/tag). see  [terra support toolkits you need](https://support.terra.bio/hc/en-us/articles/360037493971-Toolkit-All-the-tools-you-need-to-write-and-run-WDLs)


Notes:
* shell script variable. use ```${myVar}``` for all input variables. Use ```$myVar``` for local variables. Note cromwell will try and expand ```${myVar}``` even if it is in a comment

* <span style="color:red">  when you run cromwell make sure you set the user id else you will not be able to remove files created in your command task. The owner at least on mustard will be nfsnobody </span>
 

## Test a WDL file from the command line

1 how format json
```
 $ cat salmonQuantWorkflow.inputs.json | jq
{
  "quantify.salmon_paired_reads.outDir": "${}",
  "quantify.bamToFastq.runtime_preemptible": "${}",
  "quantify.salmon_paired_reads.runTimePreemptible": "${}",
  "quantify.salmon_paired_reads.memoryGb": "${}",
  "quantify.inputBam": "${this.bam_file}",
  "quantify.sampleId": "${this.sample_id}",
  "quantify.salmon_paired_reads.runTimeCpu": "${}",
  "quantify.salmon_paired_reads.dockerImg": "${}",
  "quantify.refIndexTarGz": "gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.tar.gz",
  "quantify.bamToFastq.diskSpaceGb": "${}",
  "quantify.salmon_paired_reads.diskSpaceGb": "${}",
  "quantify.bamToFastq.memoryGb": "${}"
}
```


2 validating a wdl file
cromwell and womtool are install in our conda env
```
womtool validate testBashTask.wdl
######### deprecated java -jar ../../java/bin/womtool-58.jar validate testBashTask.wdl
```

3 creating a json input template
Things to notice, Files must exist, even if you do not use them. ie just touch them

```
womtool inputs salmonPairedReadQuantTask.wdl
######### deprecated (base) [aedavids@mustard wdl]$ java -jar ../../java/bin/womtool-58.jar inputs salmonPairedRe\
adQuantTask.wdl
{
  "salmon_quant.memoryGb": "Int (optional, default = 8)",
  "salmon_quant.dockerImg": "String (optional, default = \"quay.io/biocontainers/salmon:0.14\
.1--h86b0361_1\")",
  "salmon_quant.refIndexTarGz": "File",
  "salmon_quant.sampleId": "String",
  "salmon_quant.diskSpaceGb": "Int (optional, default = 40)",
  "salmon_quant.outDir": "String",
  "salmon_quant.rightReads": "File",
  "salmon_quant.runtime_cpu": "Int (optional, default = 8)",
  "salmon_quant.leftReads": "File"
}
```

4 running cromwell
machine must have docker installed. Use mustard.

<span style="color:red">Use runCromwell.sh it sets the user id so that you can delete an files created by the task command bash script. $EUID is your user id. You can also get this by running `id -u` or just `id` </span>

```
$ ~/extraCellularRNA/bin/runCromwell.sh  --inputs salmonPairedReadQuantTask.wdl.input.json salmonPairedReadQuantTask.wdl

AEDWIP 
DEPRECATED (base) [aedavids@mustard wdl]$ ~/extraCellularRNA/bin/runCromwell.sh -jar ../../java/bin/cromwell-58.jar run --inputs salmonPairedReadQuantTask.wdl.input.json salmonPairedReadQuantTask.wdl
+ java '-Dbackend.providers.Local.config.runtime-attributes=String? docker String? docker_user="$EUID"' -jar ../../java/bin/cromwell-58.jar run --inputs salmonPairedReadQuantTask.wdl.input.json salmonPairedReadQuantTask.wdl
[2021-03-26 08:42:04,28] [info] Running with database db.url = jdbc:hsqldb:mem:9e9bff35-6e77-4505-b836-7e4a906e847f;shutdown=false;hsqldb.tx=mvcc
```

# 1vs all test
 on mustard

## step 1) create a docker image
The buildContextDir has the r script files we want to copy into the new image
 ```
cd ~/aedavids/extraCellularRNA/terra/deseq/bin

myTag="aedavids/test-1vs-all-2"
buildContextDir=../R

docker build --file ./DockerFile.1vsAll --tag $myTag $buildContextDir
 ```
 
## set the docker image tag
<span style="colore:red"> use check the docker image repository</span>
```
docker images |grep aedavids
```

set an env var
```
export IMG="aedavids/edu_ucsc_kim_lab-1vsall_1.0"
```


## step 2) test
launch the docker and open a bash shell

start the container
```
docker run --rm --detach $IMG
```

to find the name of our container
```
docker ps |grep aedavids
```

create an interactive shell on the container
```
docker exec -it nifty_lalande /bin/bash
```


to kill the container
```
docker rm -f myContainerName
```

## step 3) validate WDL

```
$ pushd ~aedavids/extraCellularRNA/java/bin

$ wget https://github.com/broadinstitute/cromwell/releases/download/74/cromwell-74.jar
$ wget https://github.com/broadinstitute/cromwell/releases/download/74/womtool-74.jar

export WDL_TOOLS=`pwd`

$ java -version
openjdk version "11.0.1" 2018-10-16 LTS
OpenJDK Runtime Environment Zulu11.2+3 (build 11.0.1+13-LTS)
OpenJDK 64-Bit Server VM Zulu11.2+3 (build 11.0.1+13-LTS, mixed mode)

$ java -jar ${WDL_TOOLS}/womtool-74.jar validate 1vsAllTask.wdl 
```

## step 4) create input and output json 
edit generated values
```
$ java -jar ${WDL_TOOLS}/womtool-74.jar inputs 1vsAllTask.wdl > 1vsAllTask.wdl.inputs.json
$ java -jar $WDL_TOOLS/womtool-74.jar outputs 1vsAllTask.wdl > 1vsAllTask.wdl.outputs.json
```



input assume we are running in a tmp sub dir
```
$ cat
cat 1vsAllTask.wdl.inputs.json 
{
  "deseq_one_vs_all.one_vs_all.countMatrix": "../../deseq/python/test/data/numReadsMatrix.csv",
  "deseq_one_vs_all.one_vs_all.design": "~ treatment",
  "deseq_one_vs_all.one_vs_all.colData": "../../deseq/python/test/data/colData.csv",
  "deseq_one_vs_all.one_vs_all.referenceLevel": "ctrl",
  "deseq_one_vs_all.one_vs_all.dockerImg": "aedavids/test-1vs-all-2",
  "deseq_one_vs_all.one_vs_all.memoryGb": "20",
  "deseq_one_vs_all.one_vs_all.diskSpaceGb": "40",
  "deseq_one_vs_all.one_vs_all.runTimeCpu": "1",
  "deseq_one_vs_all.one_vs_all.isCSV": "true"
}
```

the outputs.json file is not used by cromwell run, how ever it is useful to insure
workflows run on terra write output consistently
```
$ cat 1vsAllTask.wdl.outputs.json 
{
  "deseq_one_vs_all.one_vs_all.outfile": "deseqResult1vsAll"
}

```

## step 5) test (run cromwell) locally
<span style="color:red">THIS IS TRICKY!</span>

set isDebug='true' in 1vsAllTask.wdl.inputs.json

when you run you may see an error and stack trace like bellow. Just ignor it. your container ran. For details see [github issue](https://github.com/broadinstitute/cromwell/issues/6674)
```
[2022-03-30 15:56:09,04] [warn] BackendPreparationActor_for_6fef3552:deseq_one_vs_all.one_vs_all:-1:1 [6fef3552]: Docker lookup failed
java.lang.Exception: Unauthorized to get docker hash aedavids/edu_ucsc_kim_lab-1vsall_1.0:latest
	at cromwell.engine.workflow.WorkflowDockerLookupActor.cromwell$engine$workflow$WorkflowDockerLookupActor$$handleLookupFailure(WorkflowDockerLookupActor.scala:222)
```

You should be able to find all the output from your run under 
```
cromwell-executions/deseq_one_vs_all/randomGUID/call-one_vs_all/execution
```

runCromwell and cromwellDebug.conf set the user id so that it is easy to 
remove output files
```
cd extraCellularRNA/terra/deseq/
mkdir tmp

$ cd tmp
$ ../../../bin/runCromwell.sh \
    -Dconfig.file=../../wdl/cromwellDebug.conf \
    -jar ${WDL_TOOLS}/cromwell-74.jar run \
    --inputs ../../wdl/1vsAllTask.wdl.inputs.json \
    ../../wdl/1vsAllTask.wdl
```

check test results. File paths will have different quids
```
$ls cromwell-executions/deseq_one_vs_all/e965cf64-184a-44d4-8847-24b657823e5c/call-one_vs_all/execution/
ctrl_vs_all.results.csv  estimatedSizeFactors.csv  script.background  stderr.background
DESeqScript.out          rc                        script.submit      stdout
docker_cid               script                    stderr             stdout.background
```

## clean up images on mustard
TODO:
- docker images list
- docker images rm


# publish docker container
need to make it aviable to terra. dockerstore.org is hard to use documentation is mess and over complicated. Just use docker store, upload wdl to broad repository

- https://docs.docker.com/docker-hub/repos

on mustard

find image 1vsall image
```
$ docker images |grep aedavids | cut -f 1
aedavids/edu_ucsc_kim_lab-1vsall_1.0                                       latest                                                      6058eac32165        2 weeks ago         6.29GB
aedavids/extra_cellular_rna_2_01                                           latest                                                      c728ed04d08a        5 months ago        6.28GB
aedavids/biocontest                                                        1.0                                                         080c37b0ade1        10 months ago       3.95GB
aedavids/extra_cellular_rna_broken                                         latest                                                      b8c2a26d9690        17 months ago       5.02GB
```

create a docker repository

1. log on to
   [https://hub.docker.com/repository](https://hub.docker.com/repository)

2. press create an re repository
    * my account name is aedavids
    * new repo name is edu_ucsc_kim_lab

3. push image from mustard
   https://docs.docker.com/engine/reference/commandline/push/
  
   ```
   # image was already tag with docker user id
   $ docker login
   $ dockerUser=aedavids
   $ dockerRepo=edu_ucsc_kim_lab
   $ localTag=aedavids/edu_ucsc_kim_lab-1vsall_1.0
   $ docker image push $localTag
   the push refers to repository [docker.io/aedavids/edu_ucsc_kim_lab-1vsall_1.0]
   ca8223fbe2d3: Pushed 
   e156336e3923: Pushed 
   574ede5b3fdc: Pushed 
   21b52c92873d: Pushed 
   d3ddc5418c48: Pushed 
   f612b808c909: Pushed 
   1e138e22eaf9: Pushed 
   83160482d625: Pushed 
   c7a649503e62: Pushed 
   05e86bdff9ed: Pushed 
   5c9adbf2cdef: Pushed 
   d9608d01217a: Pushed 
   35d4453e3f63: Mounted from bioconductor/bioconductor_docker 
   c52616870693: Mounted from bioconductor/bioconductor_docker 
   8b34b2b231dd: Mounted from bioconductor/bioconductor_docker 
   8e8fdf416a80: Mounted from bioconductor/bioconductor_docker 
   84bb14749b6e: Mounted from bioconductor/bioconductor_docker 
   cf6d15ee8907: Mounted from bioconductor/bioconductor_docker 
   ae157dd031cb: Mounted from bioconductor/bioconductor_docker 
   b949cb87dde4: Mounted from bioconductor/bioconductor_docker 
   e17ef34273d7: Mounted from bioconductor/bioconductor_docker 
   5f629bbc7ac7: Mounted from bioconductor/bioconductor_docker 
   75aff22e485d: Mounted from bioconductor/bioconductor_docker 
   a75397c7fdd5: Mounted from bioconductor/bioconductor_docker 
   99384de96be4: Mounted from bioconductor/bioconductor_docker 
   7555a8182c42: Mounted from kishwars/pepper_deepvariant 
   latest: digest: sha256:330c1bcdfc532192f44f5841095936875c74bcbaf90e35e16eb54595cc68c371 size: 5997
   ```

on  https://hub.docker.com/repository/docker/aedavids/edu_ucsc_kim_lab-1vsall_1.0 you will find our image aedavids/edu_ucsc_kim_lab-1vsall_1.0

This is the matches the String dockerImg default value in 1vsAllTask.wdl

upload wdl to broads method repository. 
1. Go to terra workspace, click on 'workflows" -> find workflow -> Broad Methods Respository 
2. create new method
3 fill in form
    * namespace aedavids.ucsc.edu
    * name 1vsAllTask.wdl
    * load from file
    * export to workspace
    
