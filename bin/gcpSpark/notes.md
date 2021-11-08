
-[play list](- https://www.youtube.com/playlist?list=PLeOtIjHQdqvFtYzoFL-DCYx_Sw-5iBJQ4)


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

## second module: creating clusters, submitting jobs, ..

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

### submitting job via cloud console

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

### Submitting Job from cli
https://cloud.google.com/dataproc/docs/guides/submit-job

```
gcloud dataproc jobs submit pyspark --cluster my-frist-cluster \
    gs://aedavids-ucsc-edu-spark-tutorial2/spark_write_demo.py \
    --region us-central1 \
    --project dataprocspark-328421
    
error: the hive table already exists. Need to modify the script
write.mode("overwrite")
```

### Alternative way to submit jobs <span style="color:red"> (this is a bad idea)</span>
use standard spark-submit. Note google cloud will have an info about the job. 

ssh to master
```
spark-submit --master yarn --deploy-mode client gs://my-driver.py
```

There are other problem with this method

### deleting clusters
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
- how to size executors
- how to use preemptive vm
- local SSD
- using auto scaling





## foo
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
- https://www.youtube.com/channel/UC86J53XZybiEoYDQolBIZnw



- https://www.youtube.com/watch?v=HrGxUxcqXI8

## second module: creating clusters, submitting jobs, ..

2.1 - Create the first Dataproc Cluster | Apache Spark on Dataproc | Google Cloud Series

- ?? what does --bucket do?
I think when we did hte interactive jupyter server demo we set --bucket. this gave us access to notebooks we saved

## hacking. On the gcp dataproc create form I selected various options to see how they effect the gcloud dataproc create script

- what does --scopes do ?
Enables the cloud-platform scope for this cluster
Allow API access to all Google Cloud services in the same project.

```
gcloud dataproc clusters create my-cluster 
    --enable-component-gateway 
    --bucket aedavids-ucsc-edu-spark-tutorial2 
    --region us-central1 --zone us-central1-b 
    --master-machine-type n1-standard-2 --master-boot-disk-size 500 --num-workers 2 --worker-machine-type n1-standard-2 --worker-boot-disk-size 500 
    --image-version 2.0-ubuntu18 
    --optional-components JUPYTER 
    --scopes 'https://www.googleapis.com/auth/cloud-platform' 
    --project dataprocspark-328421
```

- we did not select 'project access' left un checked
```
gcloud dataproc clusters create my-cluster 
    --enable-component-gateway 
    --bucket aedavids-ucsc-edu-spark-tutorial2 
    --region us-central1 --zone us-central1-b 
    --master-machine-type n1-standard-2 --master-boot-disk-size 500 --num-workers 2 --worker-machine-type n1-standard-2 --worker-boot-disk-size 500 
    --image-version 2.0-ubuntu18 
    --optional-components JUPYTER 
    --project dataprocspark-328421
```
    
# enabled Personal Cluster Authentication
allows interactive we access

```
gcloud dataproc clusters create my-cluster 
    --enable-component-gateway 
    --bucket aedavids-ucsc-edu-spark-tutorial2 
    --region us-central1 --zone us-central1-b 
    --master-machine-type n1-standard-2 --master-boot-disk-size 500 
    --num-workers 2 
    --worker-machine-type n1-standard-2 --worker-boot-disk-size 500 
    --image-version 2.0-ubuntu18 
    --properties dataproc:dataproc.personal-auth.user=aedavids@ucsc.edu 
    --optional-components JUPYTER 
    --project dataprocspark-328421
```


## third module: compute storage isolation
using buckets, meta-store


## 4th performance and cost issue


## 5 moduel
- initialization
- runnign jupyter
- integration with airflow
