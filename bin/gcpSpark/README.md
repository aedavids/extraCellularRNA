# Spark related scripts

## running spark and juypter on local machine
1. use startJuypterServer.sh
2. notebook uses findspark().init to start pyspark

# creating spark clusters on GCP dataproc
I thought you could create a cluster and just stop it. How ever when you restart you can not connect to the jupyter server. you get a 500 error. It seems like the juypter server is not install permentantly

## configuring ACL (one time)

- [https://cloud.google.com/dataproc/docs/concepts/iam/iam](https://cloud.google.com/dataproc/docs/concepts/iam/iam)
Dataproc users are required to have [service account ActAs permission]( https://cloud.google.com/iam/docs/service-accounts-actas) to deploy Dataproc resources, for example, to create clusters and instantiate workflows. see [https://cloud.google.com/iam/docs/impersonating-service-accounts](https://cloud.google.com/iam/docs/impersonating-service-accounts) for details



### 1)  [understanding service accounts](https://cloud.google.com/iam/docs/service-accounts)
   * goal learn enough to get dataproce working
   * A service account is a special kind of account used by an application or a virtual machine (VM) instance, not a person. Applications use service accounts to make authorized API calls, authorized as either the service account itself, or as Google Workspace or Cloud Identity users
   * For example, a Compute Engine VM can run as a service account, and that account can be given permissions to access the resources it needs. This way the service account is the identity of the service, and the service account's permissions control which resources the service can access.
   * Service accounts do not have passwords, and cannot log in via browsers or cookies. Service accounts are associated with private/public RSA key-pairs that are used for authentication to Google. You can let other users or service accounts impersonate a service account.
   * instead of google-manged kyes. we are going to use 'user-managed keys'
   * we store both the private and public keys on our local machine. thse may be refered to as 'external keys'
   * User-managed keys can be managed by the IAM API, gcloud command-line tool, 
   * 2 types of service accounts 'user managed' and 'default service accounts'
     + 'user managed service accounts
       - name is of the form service-account-name@project-id.iam.gserviceaccount.com
       - we are responsible for CRUD
       - can have up to 100
     + 'default service accounts'
       - When you enable or use some Google Cloud services, they create user-managed service accounts that enable the service to deploy jobs that access other Google Cloud resources. These accounts are known as default service accounts.
       - names are of the form project-id@appspot.gserviceaccount.com or project-number-compute@developer.gserviceaccount.com
   * google managed service accounts
     + not listed on service accounts page in the cloud console
       - example project-number@cloudservices.gserviceaccount.com
       
### 2)  https://cloud.google.com/iam/docs/impersonating-service-accounts
* several predefined roles
  + Service Account User (roles/iam.serviceAccountUser): Allows principals to indirectly access all the resources that the service account can access.
  + Service Account Token Creator (roles/iam.serviceAccountTokenCreator):
  + Workload Identity User (roles/iam.workloadIdentityUser): Allows principals to impersonate service accounts from GKE workloads

## acl config
```
gcloud auth login

# deprecated gsutil config -a
```

## setting up gcp project specific env vars
Use the bash builtin fucntion source to run a config file. this will set environmental variable like PROJECT_ID, REGION, ... The utility scripts use these variables.

## create GCP dataproc sparc clusters

ref: [How to Create Google Cloud Dataproc Clusters with Spark + Jupyter + Python Libraries  GCP Tutorial](https://www.youtube.com/watch?v=nccCsk_MHDs)

1. create a script to install extra packages on all the notes
install.sh

2. upload install.sh to our project buck
normally you do not run upload.sh. it is part of teh createCluster.sh
upload.sh

3. run create cluster script
createCluster.sh

4. use deleteCluster.sh
destroys the cluster. No need to keep it around. we spin them up when ever we need them


# connecting to the cluster from our local machine

```
CLUSTER_NAME=my-cluster
PROJECT_ID=dataprocspark-328421

gcloud compute ssh ${CLUSTER_NAME}-m --protject=$PROJECT_ID
```


# find the URL to juypter
you can use the gcp console or 

```
CLUSTER_NAME=my-cluster
PROJECT_ID=dataprocspark-328421
REGION=us-central1

gcloud dataproc clusters describe ${CLUSTER_NAME} --region=${REGION} --project=${PROJECT_ID} | grep jupyter
```

## start pyspark notebook
ref: [Apache Spark & Jupyter on Google Cloud Dataproc Cluster  Spark + Jupyter + Dataproc](https://www.youtube.com/watch?v=5OYT2SSMGo8)

from pyspark.sql import SparkSession
spark = SparkSession.builder.getOrCreate('aedwip')


bucket = 'gs:// blah blah'

df = spark.read.csv(bucket, inferSchema=true, header=True)
