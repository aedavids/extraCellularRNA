
- [Sushil Kumar 'apache spark on dataproc Google cloud series' play list]( https://www.youtube.com/playlist?list=PLeOtIjHQdqvFtYzoFL-DCYx_Sw-5iBJQ4)


## first module intro : how to set up GCP

### setup machine for development
- install gcloud CLI

Authorization
```
gcloud auth login
Your current project is [None].  You can change this setting by running:
  $ gcloud config set project PROJECT_ID

```

 set the proeject id. I think this is easier then passing the --project
```
gcloud config set project dataprocspark-328421
Updated property [core/project].

```

I wonder if we have an account set up issue? We do not have a credentialed account? See if he deals with this in next video

```
$ gcloud auth list
Credentialed Accounts

ACTIVE: *
ACCOUNT: aedavids@ucsc.edu

To set the active account, run:
    $ gcloud config set account `ACCOUNT`
```

## 2.1 second module: creating clusters, submitting jobs, ..

create a cluster using the gcp consol  data proc create form. very trivial cluster. no auto scaling. did not enable component gateway, or select any options like jupyter

manage seucrity: "allow api access to all google cloud services in the same project"

google-managed encrption keys
```
gcloud dataproc clusters create my-frist-cluster \
    --region us-central1 --zone us-central1-a \
    --master-machine-type n1-standard-2 --master-boot-disk-size 500 \
    --num-workers 2 --worker-machine-type n1-standard-2 \
    --worker-boot-disk-size 500 \
    --image-version 2.0-ubuntu18 \
    --scopes 'https://www.googleapis.com/auth/cloud-platform' \
    --project dataprocspark-328421
```

cluster details -> web interfacess provides script for creating an ssh tunnel

### 2.2 submitting job via cloud console

create a bucket to store results
```
gsutil mb gs://my-globally-uniqu-name
```

download sample driver python program https://gist.github.com/kaysush/f842155dfbee07701adea236b491fbab

copy to a bucket. We should be able to use file:// how ever I got a permission problem

```
$ gsutil cp ./spark_write_demo.py   gs://${BUCKET}/spark_write_demo.py

$ gsutil ls gs://${BUCKET}
gs://aedavids-ucsc-edu-spark-tutorial2/install.sh
gs://aedavids-ucsc-edu-spark-tutorial2/rose.txt
gs://aedavids-ucsc-edu-spark-tutorial2/spark_write_demo.py
```

fill in submit form. Notice property spark.submit.deploymode. We wnat to use client because we have a deicated master. Leave works free to do as much work as possible

```
DO NOT USE REST
POST /v1/projects/dataprocspark-328421/regions/us-central1/jobs:submit/
{
  "projectId": "dataprocspark-328421",
  "job": {
    "placement": {
      "clusterName": "my-frist-cluster"
    },
    "statusHistory": [],
    "reference": {
      "jobId": "my-first-job-2",
      "projectId": "dataprocspark-328421"
    },
    "pysparkJob": {
      "mainPythonFileUri": "gs://aedavids-ucsc-edu-spark-tutorial2/spark_write_demo.py",
      "properties": {
        "spark.submit.deploymode": "client"
      }
    }
  }
}
```

check the results. we wrote to a hive table

set jobType to Hive

set quer source to query text

set query text to '
```
DO NOT USE REST!
POST /v1/projects/dataprocspark-328421/regions/us-central1/jobs:submit/
{
  "projectId": "dataprocspark-328421",
  "job": {
    "placement": {
      "clusterName": "my-frist-cluster"
    },
    "statusHistory": [],
    "reference": {
      "jobId": "my-hive-job",
      "projectId": "dataprocspark-328421"
    },
    "hiveJob": {
      "properties": {},
      "scriptVariables": {},
      "queryList": {
        "queries": [
          "select * from random_numbers"
        ]
      }
    }
  }
}
```

submit. results from output window

```
+--------------------+
| random_numbers.id  |
+--------------------+
| 50                 |
| 51                 |
| 52                 |
| 53                 |
| 54                 |
| 55                 |
| 56                 |
```

data is on cluster. If we delete cluster we will loose data. will deal with this latter

### 2.3 Submitting Job from cli
https://cloud.google.com/dataproc/docs/guides/submit-job

```
gcloud dataproc jobs submit pyspark --cluster my-frist-cluster \
    gs://aedavids-ucsc-edu-spark-tutorial2/spark_write_demo.py \
    --region us-central1 \
    --project dataprocspark-328421
    
error: the hive table already exists. Need to modify the script
write.mode("overwrite")
```

### 2.4 Alternative way to submit jobs <span style="color:red"> (this is a bad idea)</span>
use standard spark-submit. Note google cloud will have an info about the job. 

ssh to master
```
spark-submit --master yarn --deploy-mode client gs://my-driver.py
```

There are other problem with this method

### 2.5 deleting clusters
we have a script for this

In our example we write to a hive table. By default they are stored on cluster. If we delete cluster we loose data and hive meta data

we also loose everything writtent to HDFS (hadoop file system)

<span style="color:red">If we stop the cluster we will loose all of our data</span>

## Module 3 storage and compute isolation

###write to google cloud storage bucket
we want to be able to delete cluster with out loosing data. running the cluster is expensive

move storage to google cloud storage (gcs). Just change write() pass bucket url, All spark gcp instannce come preconfigured to read/write to gcs

### externalize HIVE storage
This is long viedo we do not use hive ref: https://www.youtube.com/watch?v=dfq5RuQ32HI&list=PLeOtIjHQdqvFtYzoFL-DCYx_Sw-5iBJQ4&index=10


google cloud console -> SQL

create mySql instance
- no password
- put in same region as cluster

## Module 4: performance

### 4.1 resource and performance optiomization

### 4.2  use trainient per job clusterms
- better fit, more cost effective

### 4.3 sizing the executor
- YARN by design only able to use 80% of memory
- tune number of work for level of parallization (ie number of executors)
- example 1
  + cluster with 10 works, (each mchine 16 vcores and 64 GB memory )
  + total resource = 160 vcores and 640 GB
  + 80% of memory == 500 GB of memory
  + each worker == 50 GB and 16 vcores
  + 50 GB / 16 vCores == 3.125 GB / vCore. To make math easy examples assumes 3 not 3.125
  + if we creat executors with 4 vcores and 12 GB we can configure executors in 2 ways
    * 1) 11 GB for executor mem and 4 vCores with 1 GB for over head 
    * 2) 10 GB for executor mem and 4 vcores with 2 GB for over head
  + we packed 4 executor in a work node
    * 4 executors * 10 works == 40 executors
    * 160 total vcores
    * can run upto 160 parallel tasks
  + <span style="color:red">if max parralization of job is only 100 our cluster can not use all the resource reduce the number of workers</span>
- example 2
  + high memory machine (16 vcores and 128 GB memroy for each machine) 
  + 10 workers == total reasource 160 vcore and 1.2 TB memroy
  + 80% of memory is about 1 TB GB for all 10 workers
  + total vcores = 160
  + each worker == 100 GB and 16 vcores
  + 6 GB / 1 vcore
  + if we configure 24 GB and 4 vcore per executor 2 choices
    * 22 GB and 4 core with 2 gb for over head or 20 gb adn 4 vcores with 4gb for over head
    
- example 3 counter example
  + 16 vcores per work
  + create exector with 5 workers 
    * we can only use at most 15 vcores 
    * <span style="color:red">we waste 1 vcore per worker! if we have 10 works in cluster we waste 10 vcores</span>

- example 4 counter example
  + 50 gb per worker
  + if we use more than 12 gb for each executor we will not be able to fit 4 executor, some vcodes will not be used 
  
### 4.4 preemptibile secondary workers
- what are preemptible VMs
  + like regular VM but can be taken away at any time
  + will be taken away after 24 hrs of use
  + much cheaper. 
  
- compute only works
- same machine type as primary workers
- do not host HDFS. 
  + not an issue for us we only use gcp storage
- spark is resilient, if preemptition event occures spark will reschedule task
- caveats
  + do not use more than 50% of your total works as preemptible
  + may experience high number of task failures due to preemption
  + need to make cluste more tolerable to failures
    * set parameters
      ^ yarn:yarn.resouremanager.am.max-attempts=10
      ^ spark:spark.task.maxFailures=10
      ^ spark:spark.stage.maxConsecurtiveAttemps=10

- tip
  - get a baseline wiht primary works only
  - add preemptible secondary works on top
  - finish jobs well with in SLA with total lower cost
  
### 4.5 local SSD for resecue
- by default works on have boot disk
- spark uses to write shuffle files
- add local SSD they are fast and cheap
- if you are using secondary workers make sure you add sdd to them

### 4.6 using auto scaling
- based on amount of pending memory vs aviable memory
- configure scaleUp or scaleDown formula value between 0 and 1 
- frequency control (how often dataproc check scaling policy)
- scaleUpMinWorkderFraction and scaleDownMinWorkFraction
- select autoscaling policies -> create -> spark with dynamic allocation
  * this policy is good 
  * policies are region specific
  * DO NOT SCALE PRIMARY works that use HDFS !!!
  
- use autoscaling when storge is externalized. 
  + i.e. not HDFS
  + many jobs running on a single cluster
  + single job wit variable resource requirements
  
- never use autoscaling
  - if you use HDFS or spark streaming (2018)
  - for implementing idle cluster
    - just use tranisent clusters
    
- if you use auto scaling make your cluster more fault tolerent
    * set parameters
      ^ yarn:yarn.resouremanager.am.max-attempts=10
      ^ spark:spark.task.maxFailures=10
      ^ spark:spark.stage.maxConsecutiveAttempts=10




## 5 moduel TODO watch
- initialization
- runnign jupyter
- integration with airflow
  - <span style="color:red">airflow is an apache workflow that runs across different platforms. it spark could be a task </span>
  
### 5.1 Juypter on Dataproc
- create cluster
- select 'enable component gateway'
- select Jupyter notebook component
- by default juypter notebook will be stored in staging bucket
  - to change where notebook is stored add property
    dataproc jupyter.notebook.gcs.dir myBucket/notebooks
- finish cluster create
- go to cluster page
- select tab web interface
- you will see link for jupyter and juypter lab
- go to vm instance tab
  - get public ip address of master
  
#### finding URL from CLI
```
gcloud dataproc clusters describe ${CLUSTER_NAME} --region=${REGION} --project=${PROJECT_ID} | grep jupyter
clusterName: jupyter-test1
      Jupyter: https://i32ejq3punf7toyrbo3hm7qe3u-dot-us-central1.dataproc.googleusercontent.com/jupyter/
      JupyterLab: https://i32ejq3punf7toyrbo3hm7qe3u-dot-us-central1.dataproc.googleusercontent.com/jupyter/lab/
    - jupyter-test1-m
      dataproc:jupyter.notebook.gcs.dir: anvil-gtex-v8-hg38-edu-ucsc-kim-lab-spark/notebooks
      hdfs:dfs.namenode.lifeline.rpc-address: jupyter-test1-m:8050
      hdfs:dfs.namenode.servicerpc-address: jupyter-test1-m:8051
    - jupyter-test1-w-0
    - jupyter-test1-w-1
  goog-dataproc-cluster-name: jupyter-test1
``

