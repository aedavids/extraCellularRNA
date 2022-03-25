# ref
- [google cloude forum dataproc label](https://www.googlecloudcommunity.com/gc/forums/filteredbylabelpage/board-id/cloud-data-analytics/label-name/dataproc)

- [apache spark faq + email list info](https://spark.apache.org/faq.html)

- [apache spark community stackOverflow, mailing list](https://spark.apache.org/community.html)
- [apache spark config](https://spark.apache.org/docs/latest/configuration.html#spark-configuration)

  * view all properties, executorion behavior, logging. env vars, ...
- [building production PySpark jobs](https://medium.com/@lubna_22592/building-production-pyspark-jobs-5480d03fd71e)

- [spark join strategies](https://towardsdatascience.com/strategies-of-spark-join-c0e7b4572bcf)

- [gcp cluster web interfaces (YARN ...)](https://cloud.google.com/dataproc/docs/concepts/accessing/cluster-web-interfaces)
  * 3 ways
    1. use dataproc component gateway 
       - this is how we run juypter notebooks
    2. cloud shell. 
       - use gcp cluster console 'web preview' 
    3. gcloud
       - gcloud compute ssh to create ssh tunnel
       ```
       gcloud compute ssh ${HOSTNAME} \
           --project=${PROJECT} --zone=${ZONE}  -- \
           -D ${PORT} -N
        ```
        
- [monitoring spark web interfaces, and after the fact history](https://spark.apache.org/docs/latest/monitoring.html)

- [gcp dataproc doc](https://cloud.google.com/dataproc/docs)
- [gcp doc](https://cloud.google.com/docs)
- spark performancs optimation
  * [skew # 1](https://medium.com/road-to-data-engineering/spark-performance-optimization-series-1-skew-2762a0f288c)
  * [spill # 2](https://medium.com/road-to-data-engineering/spark-performance-optimization-series-2-spill-685126e9d21f)
  * [shuffle](https://medium.com/road-to-data-engineering/spark-performance-optimization-series-3-shuffle-104738a83a9e)
  
# configure logging
seems like by default our cluster logs get deleted when the cluster deactivate

* where are the driver and worker logs?
* [dataproc quides/logging](https://cloud.google.com/dataproc/docs/guides/logging)
   + [cloud logging priceing summary](https://cloud.google.com/stackdriver/pricing)
   + free per month; first 50 GiB
   + $0.01/Gib for logs retained more than 30 days
   + logs retained for default retention period do not incur a storage cost
   + log bucket _Default ; 30 day retention period
 * [Dataproc job logs in Logging](https://cloud.google.com/dataproc/docs/guides/logging?authuser=1#job_logs_in)
   + enable job dirver logs in cloud logging
     ```
     dataproc:dataproc.logging.stackdriver.job.driver.enable=true
     ```
   
   + required properties set by default when cluster is created
     ```
     dataproc:dataproc.logging.stackdriver.enable=true
     dataproc:jobs.file-backed-output.enable=true
     ```
   
   + Enabling YARN container logs in Cloud Logging
     ```
     dataproc:dataproc.logging.stackdriver.job.yarn.container.enable=true
     ```
   
   + required property set by default when a cluster is created;
     ```
     dataproc:dataproc.logging.stackdriver.enable=true
     ```
   
 * [accessing logs in cloud logging](https://cloud.google.com/dataproc/docs/guides/logging?authuser=1#accessing_job_logs_in)
   + Dataproc Job driver and YARN container logs are listed under are listed under the Cloud Dataproc Job resource.
   + can also access using 'log explore' service
   + can also use gcloud to read to local machine



# where does it run out of memory?
seems to blow up with  same amount of elapsed time

parts of library normilization require entire data set to be in memory

## Create a small crash test on GCP
- Make it faster and easier to debug and test
- Use 1 master 2 workers
- try 100 quant files

## where does OOM origininate?
- driver ?
- coloces ?
- executor runs out of memory? several parts to executor memory

## TODO
- are we using logging in python correctly?

- is spark using the extra local SSD?

- create debugLaunchCluster
  - do we have to do something to keep history and logs
  - open ssh tunnel for yarn monitoring
    - do we spill?
    
- DONE change logger info -> warn

- DONE run driver on worker
  * did not make a difference

- use more local machine

- where do worker logs go?

- HDFS is several TB. utilization is zero???


- ? Sizing and configuraiton ?
  * where are we crashing?
  * are we able to create fewer executors with more memory?

- run wiht one m1 worker

- do we nee do anything to make it easier to work with log files when we nave multiple worker
  https://medium.com/@lubna_22592/building-production-pyspark-jobs-5480d03fd71e
  
- display query plan
  - does it change as it grows?
  - BroadcastHashJoin exhaust memory 
    - https://towardsdatascience.com/strategies-of-spark-join-c0e7b4572bcf
    
    https://github.com/apache/spark/blob/master/sql/core/src/main/scala/org/apache/spark/sql/execution/SparkStrategies.scala#L111
    
- can we avoid join
  we can select a col and append it
 bad  https://stackoverflow.com/questions/42853778/add-a-column-from-another-dataframe
 !!!! this is the answer https://sparkbyexamples.com/spark/spark-add-new-column-to-dataframe/ 
  
- if my rows are ordered will the part files be ordered?

- what if we do not collece before writting. Lots of parts file use shell script to concatinte them

# partition quant file design
  * split quant files into x parts.
    + see unix split command
    + the first split for each quant file must have the same row names
  * for each split create a raw count dataframe
  * use spark union to combine the raw counts into master raw counts
  * run library compostition normalization ie calcualte scaling factors
  
  * https://stackoverflow.com/a/37849488/4586180
  * https://sparkbyexamples.com/pyspark/pyspark-partitionby-example/

# Write a C program (low memory file based map reduce)
It is not clear where the OOM execption comes form? I would expect spark would spill partistions durring join/map reduce when memory gets tight how ever we never see HDFS usage change. It is always zero.



## map stage
```
For quantFileName in listOfQuantFiles
    for each line
        select transcriptName, numReads
        write( fileName=quantFileName/nametranscriptName, value = numReads )
```

## reduce stage
```
for transcripName in listOfTranscripts
    for each quant file
        numReads = read( file=quantFile/transcriptName )
        str = numReads + ","
        write(file=rawCounts, str)
    write(file=rawCounts, "\n" )
```

This is a lot of file I/O how ever if we use SSD might not be to bad

## running concurrently
1. You can split the listOfQuantFiles into separate map stage batches
2. you can split the listOfTranscripts into seperate reduce stage batches.
   - make sure each batch writes a separate file
   - last step is to cat/spark union the files into a single file

## performance considerations
- We need to avoid writting millions of tiny files to google cloud storage buckets
- map stage batches
  - local quant files (copy to local file system
  - write to local file system
  - create zip file
  - move zip file to google cloude storage buckets
- reduce stage
  - create batch specific zip file containing <span style="color:red">AEWIP</span>
  - select TODO


