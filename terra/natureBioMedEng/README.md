# natureBioMedEng
```
aedavids@ucsc.edu
2/28/23
```

**TODO: see extraCellularRNA/deconvolutionAnalysis/bin/pipeline.sh for example of how to use the cromwell run  --options argment to copy output files to a known location**

/private/groups/kimlab/natureBioMedEng

Use this directory to keep track of work we did for the manuscript resubmission deadline (march 8th) of  our nature biomedical engineering manuscript

Daniel: "we want to"

* 1  DESeq2 using all TCGA pancreatic cancer vs. all GTEx pancreas RNA-seq samples
    + i.e. PAAD vs pancreas

* 2 Whole blood vs GTEx. do not include TCGA

ref:
- extraCellularRNA/terra/wdl/README.md
- extraCellularRNA/terra/deseq/R/README.md
- https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#model-matrix-not-full-rankzy


***
# GTEx WholeBlood vs GTEx

our wdl runner, cromwell is written in Java. The only thing we use from the conda extraCellularRNA enviroment is openjdk. Mustard has an old version intalled

```
aedavids@mustard $ which java
/usr/bin/java

aedavids@mustard $ java -version
openjdk version "1.8.0_171"
OpenJDK Runtime Environment (build 1.8.0_171-b10)
OpenJDK 64-Bit Server VM (build 25.171-b10, mixed mode)

aedavids@mustard $ conda activate extraCellularRNA
(extraCellularRNA) which java
~/miniconda3/envs/extraCellularRNA/bin/java

(extraCellularRNA) java --version
openjdk 11.0.1 2018-10-16 LTS
OpenJDK Runtime Environment Zulu11.2+3 (build 11.0.1+13-LTS)
OpenJDK 64-Bit Server VM Zulu11.2+3 (build 11.0.1+13-LTS, mixed mode)
```

```
conda activate extraCellularRNA

(extraCellularRNA) pwd
/private/home/aedavids/extraCellularRNA/terra/natureBioMedEng

$ mkdir GTExWholeBlood.vs.GTEx
cd  GTExWholeBlood.vs.GTEx

$ cp ../test/test.1vsAllTask.input.json GTExWhole_Blood.vs.GTEx.1vsAllTask.input.json

$ export WDL_TOOLS=/private/home/aedavids/extraCellularRNA/java/bin
```

```
(extraCellularRNA) cat GTExWhole_Blood.vs.GTEx.1vsAllTask.input.json 
{
    "deseq_one_vs_all.one_vs_all.memoryGb": "64",
    "deseq_one_vs_all.one_vs_all.runTimeCpu": "2",
    "deseq_one_vs_all.one_vs_all.colData": "/private/groups/kimlab/GTEx/GTExTrainColData.csv", 
    "deseq_one_vs_all.one_vs_all.isCSV": "true",
    "deseq_one_vs_all.one_vs_all.isDebug": "true",
    "deseq_one_vs_all.one_vs_all.referenceLevel": "Whole_Blood",
    "deseq_one_vs_all.one_vs_all.diskSpaceGb": "80",
    "deseq_one_vs_all.one_vs_all.design": "~  sex + tissue_id",
    "deseq_one_vs_all.one_vs_all.dockerImg": "aedavids/edu_ucsc_kim_lab-1vsall_1.0",
    "deseq_one_vs_all.one_vs_all.runTimePreemptible": "1",
    "deseq_one_vs_all.one_vs_all.countMatrix": "/private/groups/kimlab/GTEx/GTExTrainGroupByGenesCountMatrix.csv"
}
```

## Run

```
touch run.sh
chmod u+x run.sh
```

Create an easy of use script
```
$cat run.sh
cat run.sh 
#!/bin/sh
# aedavids@ucsc.edu
# 3/1/23

# script makes it easy to run batch job

set -x # turn debug trace on

../../../bin/runCromwell.sh \
    -Dconfig.file=../../wdl/cromwellDebug.conf \
        -jar ${WDL_TOOLS}/cromwell-85.jar \
        run \
        --inputs ./GTExWhole_Blood.vs.GTEx.1vsAllTask.input.json \
        ../../wdl/1vsAllTask.wdl
```

use setsid to run docker in a new session with out a terminal. This will allow us
to log out with out killing the job

```
setsid sh -c 'set -x; run.sh' > run.sh.out 2>&1 &
```

monitoring batch jobs

use ps
```
userId=aedavids
$ ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep $userId
```

use pstree
```
$ pstree $userId
```


killing all the children
<span style="color:red">the syntax is strange notice the ‘-’</span>
```
kill -TERM -[pgid]
```

### debugging
debug if your job fails look under cromwell-exectution/ for DESeqScript.out, stdout and stderr . chances are your design is wrong

If you job completed sucesfully you will the total run time 

```
aedavids@mustard $ tail output/DESeqScript.out 
3       (CA)n 1437.09144        1.76210 0.0444037         0         0
4     (CAAC)n    5.99613        5.79705 0.1040761         0         0
5    (CATCT)n   69.19080       -5.11043 0.1015376         0         0
6     (CCCT)n    4.60658        2.45073 0.0548266         0         0
debug callign write.table
wrote file:  Whole_Blood_vs_all.results.lfcShrink 

*********DEBUG end saveResults()DESeqScript.R completed sucessfully

run timeTime difference of 1.026396 days
```

## copy output and remove cromwell log file and localized data
cromwell create a lot of log files. It also copied the colData and groupbyGene data file under the cromwell-executions/. We want to copy the output files and delete everything else

list all the files in all the sub dir
```
$ find cromwell-executions -type f 
```

copy the results, Note the paths will vary from run to run

```
$ natureBioMedEng/GTExWhole_Blood.vs.GTEx
$ mkdir output

$ cp cromwell-executions/deseq_one_vs_all/2479b704-2c04-42c8-8fc0-d924e3f79c87/call-one_vs_all/execution/estimatedSizeFactors.csv output

$ cp cromwell-executions/deseq_one_vs_all/2479b704-2c04-42c8-8fc0-d924e3f79c87/call-one_vs_all/execution/Whole_Blood_vs_all.results output

$ cp cromwell-executions/deseq_one_vs_all/2479b704-2c04-42c8-8fc0-d924e3f79c87/call-one_vs_all/execution/Whole_Blood_vs_all.results.lfcShrink output

# you do not need to copy the following log file
# is document exactly what happed
$ cp cromwell-executions/deseq_one_vs_all/2479b704-2c04-42c8-8fc0-d924e3f79c87/call-one_vs_all/execution/DESeqScript.out output
```

```
aedavids@mustard $ ll output/
total 14M
-rw-r--r-- 1 aedavids prismuser  13K Mar  2 20:44 DESeqScript.out
-rw-r--r-- 1 aedavids prismuser 6.0M Mar  2 20:41 Whole_Blood_vs_all.results.lfcShrink
-rw-r--r-- 1 aedavids prismuser 7.0M Mar  2 20:40 Whole_Blood_vs_all.results
-rw-r--r-- 1 aedavids nfsnobody 177K Mar  1 19:32 estimatedSizeFactors.csv
```

clean up
```
$ rm -rf cromwell-executions cromwell-workflow-logs run.sh.out 
```

***
# PAAD.vs.pancrease

use natureBioMedEng/PAAD.vs.pancrease/createDatasetForPAADvsPancreas1vsAllAnalysis.ipynb create colData and groubyGeneCountMatrix. (there is an html version)

input json
```
aedavids@mustard $ cat PAAD.vs.pancrease.1vsAllTask.input.json 
{
    "deseq_one_vs_all.one_vs_all.memoryGb": "64",
    "deseq_one_vs_all.one_vs_all.runTimeCpu": "2",
    "deseq_one_vs_all.one_vs_all.colData": "/private/groups/kimlab/natureBioMedEng/PAAD.vs.pancrease/pancreasePAADColData\
.csv",
    "deseq_one_vs_all.one_vs_all.isCSV": "true",
    "deseq_one_vs_all.one_vs_all.isDebug": "true",
    "deseq_one_vs_all.one_vs_all.referenceLevel": "Pancreas",
    "deseq_one_vs_all.one_vs_all.diskSpaceGb": "80",
    "deseq_one_vs_all.one_vs_all.design": "~  Gender + Cohort",
    "deseq_one_vs_all.one_vs_all.dockerImg": "aedavids/edu_ucsc_kim_lab-1vsall_1.0",
    "deseq_one_vs_all.one_vs_all.runTimePreemptible": "1",
    "deseq_one_vs_all.one_vs_all.countMatrix": "/private/groups/kimlab/natureBioMedEng/PAAD.vs.pancrease/pancreasePAADGro\
upByGenseCountMatrixData.csv"
}
```

create easy of use script
```
$pwd
terra/natureBioMedEng/PAAD.vs.pancrease
$ cat run.sh 
#!/bin/sh
# aedavids@ucsc.edu
# 3/1/23

# script makes it easy to run batch job

set -x # turn debug trace on

../../../bin/runCromwell.sh \
    -Dconfig.file=../../wdl/cromwellDebug.conf \
	-jar ${WDL_TOOLS}/cromwell-85.jar \
	run \
	--inputs ./PAAD.vs.pancrease.1vsAllTask.input.json \
	../../wdl/1vsAllTask.wdl
```


## Run
```
$ conda activate extraCellularRNA

$ export WDL_TOOLS=/private/home/aedavids/extraCellularRNA/java/bin

$ setsid sh -c 'set -x; run.sh' > run.sh.out 2>&1 &
```

cromwell-executions/deseq_one_vs_all/6c1f29a9-d677-4b97-9ece-3ebf8dc04e92/call-one_vs_all/execution/stdout

converting counts to integer mode
Error in checkFullRank(modelMatrix) : 
  the model matrix is not full rank, so the model cannot be fit as specified.
  One or more variables or interaction terms in the design formula are linear
  combinations of the others and must be removed.
