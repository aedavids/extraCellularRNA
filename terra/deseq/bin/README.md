# TODO docker
publish image
- https://docs.docker.com/docker-hub/
- https://hub.docker.com/repositories
- https://docs.dockstore.org/en/stable/getting-started/getting-started.html
- https://docs.dockstore.org/en/stable/getting-started/register-on-dockstore.html



# AEDWIP how big a cluster do we need?

our salmon TE aware ref index has  5387495 unique transcripts. Once I group the transcripts by gene id we get 74777 unique genes. The GTEx train set has 10411 samples

- GTExTrainNumReadsMatrix.tsv 
```
ll GTExTrainNumReadsMatrix.tsv
-rw-r--r-- 1 aedavids prismuser 213G Jan 19 20:19 GTExTrainNumReadsMatrix.tsv

# number of rows
$ wc -l GTExTrainNumReadsMatrix.tsv
5387496 GTExTrainNumReadsMatrix.tsv

# number of numeric columns
$ head -n 1 GTExTrainNumReadsMatrix.tsv | awk -F$'\t' '{print NF-1;}'
10408


#numeric memory requirement
# https://medium.com/@robertopreste/tabular-data-memory-requirements-9881d2bf747a
>>> 5387496 * 10408 * 8 / 2**20/1024
417.77 GB

# memory required for transcript names
$ wc -l transcriptNames.txt 
5387496 transcriptNames.txt
$ ll transcriptNames.txt
-rw-r--r-- 1 aedavids prismuser 470M Jan 19 20:02 transcriptNames.txt

# memory requirements for sample names
$ head -n 1 GTExTrainNumReadsMatrix.tsv | wc -c
262607
```


- countsGroupedByGene/part*.csv
```
# numeric memory only
>>> 74777 * 1040 * 8 / 2**20/1024
5.7941734790802 GB
```


# List of files

- copyQuantFilesFromTerra.sh
  used to create count matrix. Copies number of Reads files from gcp to local machine.
  
- createBigDataDeseqZip.sh
  creates a python archive. Used by submit script to ensure jobs run latest versions

- createGTExTrainingCluster.sh
  create a GCP dataproc cluster running apache spark. Edit script to change cluster size.
  takes a configuration file as an arguments. Ex. data/GTEx-training.config

- createSalmonNumReadsMatrix.sh
  easy of use wrapper around transpose.py. Failed attempt to create count matrix by first combining the numReads columns into a transposed form of the count matrix using sed. row id are the salmon transcript refernece index. traspose.py reads was unable to use pandas to read file and transpose into correct deseq format. ie. columns are samples
  
  see createSalmonNumREadsTransposeMatrix.sh
  
- createSalmonNumReadsTransposeMatrix.sh
    see createSalmonNumReadsMatrix.sh

- createSparkQuantFileFromTerraQuantFiles.sh
  moves quant files from Terra to a native GCP project that can run dataproc apache spark clusters

- dateStamp.sh
  creates a time stamp. Usefull for naming output files

- diagnoseCluster.sh
  runs a gcp dataproce utility script. collection information useful for debugging cluster problems

- downloadNumReads.sh
  easy of use script for moving files form gcp bucket to local computer

- parseSalmonPartitionReadsCLI.sh,  parseSalmonReadsSelectCountsCLI.sh, parseSalmonReadsSplitCLI.sh
  wrappers around python programs

- runGTExTrain.sh
  create cluster and run job

- startSparkHistoryServe.sh
  gcp data proc does not run the spark web ui. This makes debugging difficult. They do support running the sparkHistory servers. See description bellow
  
- submitSparExtimateScalingFactors.sh <span style="color:red">use to group counts by genes</span>
  - input is a tsv matrix. Each column is the numReads col result from salmon quantification. calculates the estimated scaling factors and performance group by geneId

-  preprocessData.py. does not work on big data sets. See 'How GTExSalmonReadsMatrix.csv was created' bellow



# how to run test from cli
these are not unit test

```
cd extraCellularRNA/terra/deseq/bin
```

update the python zip file
```
createBigDataDeseqZip.sh
```

set env vars
```
pushd ../spark-3.1.2-bin-hadoop3.2/
export SPARK_HOME=`pwd`
popd
pushd ../python
export PYTHONPATH="${PYTHONPATH}:${SPARK_HOME}/python:./bigDataDeseq"
```

start conda env. 
```
sd
cap
cd $d
```

run
```
python bigDataDeseq/preprocessData.py --help

python bigDataDeseq/preprocessData.py --quantFilesCSV=test/data/mockQuantFiles.csv --mappingCSV=test/data/mockTxId2GeneId.csv --outputDir=tmp
```

# How to view Spark History locally
gcp does not provide access to the spark web UI. How ever if a completes, crashes or is stop (<span style="color:red">delete removes history and log files</span>) you can down load the history file and display them locally

### step 1 downloads
1. log on to dataproc cluster detail web page
2. select web interfaces
3. select spark history server
4. download the from history server

### step 2 start history server on local machine
1. set variable to location of download on local machine
```
historyDir=pathToDownload
export SPARK_HISTORY_OPTS="-Dspark.history.fs.logDirectory=file:${historyDir}";
```

2. start the server
```
$SPARK_HOME/sbin/start-history-server.sh
open http://localhost:18080
```

### step 3 stop the history server
```
$SPARK_HOME/sbin/stop-history-server.sh
```

## test parseSalmonReadsSplitCLI.py
test on locally
```
d=~/googleUCSC/kimLab/extraCellularRNA/terra/deseq/bin/data/spark6SampleTest.out
python bigDataDeseq/parseSalmonReadsSelectCountsCLI.py -q test/data/mockQuantFiles.csv -o  $d

# create a batch file
batchFile=~/googleUCSC/kimLab/extraCellularRNA/terra/deseq/python/test/data/parseSalmonReadsSelectCountsCLI-batch-0

pushd $d

ls  parseSalmonReadsSelectCountsCLI.out/*/p*.csv > t

echo "name,source" > $batchFile
for f in `cat t`;
do
    sampleName=`echo $f | cut -d / -f 2`
    echo "${sampleName},${d}/${f}" >> ${batchFile}
done


popd

python bigDataDeseq/parseSalmonReadsSplitCLI.py -n 2 -o $d -q ${batchFile}
```


# How GTExSalmonReadsMatrix.csv was create
We where unable to use Apache Spark to select the numReads column from each quant.sf.gz file and then join all the columns together. See bellow for more details. Our work around was to use spark to select the numReads columns from each quant.sf file and write them to a gcp bucket. Spark made it easy to read a gz file from a bucket. The file name is the sample name. The file has a single column. The name of the column is the sample name

Next we copyied the files to mustard:/scratch/aedavids. this machine has something like 1 Tb of Ram. we used the unix paste command to create the final count matrix

## detailed instructions

1. use spark to select the numReads columns from each quant.sf file
see  See terra/deseq/bin/parseSalmonReadsSelectCountsCLI.sh

2. copy the numReads files to mustard /scratch/aedavids
mustard is a machine with a lot of memory

    a. set up some basic configuration variables
        ```
        $ pwd
        /private/home/aedavids/extraCellularRNA/terra/deseq/bin/data
        (extraCellularRNA) aedavids@mustard $ source GTEx-training.config 
        ++ export BUCKET=anvil-gtex-v8-hg38-edu-ucsc-kim-lab-spark
        ++ BUCKET=anvil-gtex-v8-hg38-edu-ucsc-kim-lab-spark
        ++ export REGION=us-central1
        ++ REGION=us-central1
        ++ export PROJECT_ID=dataprocspark-328421
        ++ PROJECT_ID=dataprocspark-328421
        ++ export ZONE=us-central1-c
        ++ ZONE=us-central1-c
        ```

    b. create list of files we need to copy
        ```
        $ cd /scratch/aedavids/
        $ gsutil ls gs://${BUCKET}/quant/sparkGTEx-train.out/parseSalmonReadsSelectCountsCLI.out/onlyNumReads/ > numReadsFileList.txt

        # expected number of files in the GTEx training set
        $ wc -l numReadsFileList.txt
        10411 numReadsFileList.txt
        ```
        
    c. split numReadsFilesList into three patches

        ```
        $ split numReadsFileList.txt --numeric-suffix --lines=3000 numReadsFileListBatch.txt.

        $ ls -1 numReadsFileList*
        numReadsFileListBatch.txt.00
        numReadsFileListBatch.txt.01
        numReadsFileListBatch.txt.02
        numReadsFileListBatch.txt.03
        numReadsFileList.txt


        # sanity check
        $ wc -l numReadsFileList*
        2594 numReadsFileListBatch.txt.00
        2596 numReadsFileListBatch.txt.01
        2602 numReadsFileListBatch.txt.02
        2619 numReadsFileListBatch.txt.03
        10411 numReadsFileList.txt
        ```
 
     d. run downloads use gsutil rsync to recover from errors
        We want to make sure jobs contine to run even if we log out


        use setsid 
    
        This will let your batch job continue running after you log out.
        setsid will make it easy to kill your bulk download process and any child process. 
        See extraCellularRNA/deconvolutionAnalysis/bin/best25GTEx_TCGA.sh and  extraCellularRNA/bin/exRNADownload/batchDownloadExRNA.orgData.sh for example. There are some trick for passing arguments 
        to the child processes
        ```
        $ setsid sh -c 'set -x; myBatchJob.sh ' > $scriptLog 2>&1 &
        ```
        
    e. monitoring your batch jobs using log files
        if redirected stdout and stderr as in the above you can use the tail -f command

        use ps to montor your jobs
        ```
        userId=aedavids
        $ ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep $userId
        ```

        use pstree
        ```
        (base) [aedavids@plaza bin]$ pstree $uid
        ```

        killing all the childred. <span style="color:red">the syntax is strange notice the '-' </span>
        ```
        kill -TERM -[pgid]
        ```


        ```
        setsid sh -c 'time /private/home/aedavids/extraCellularRNA/terra/deseq/bin/downloadNumReads.sh numReadsFileListBatch.txt.00 sparkGTEx-train.out/parseSalmonReadsSelectCountsCLI.out/onlyNumReads ' > downloadNumReads.sh.batch.00.out 2>&1 &

        setsid sh -c 'time /private/home/aedavids/extraCellularRNA/terra/deseq/bin/downloadNumReads.sh numReadsFileListBatch.txt.01 sparkGTEx-train.out/parseSalmonReadsSelectCountsCLI.out/onlyNumReads ' > downloadNumReads.sh.batch.01.out 2>&1 &

        setsid sh -c 'time /private/home/aedavids/extraCellularRNA/terra/deseq/bin/downloadNumReads.sh numReadsFileListBatch.txt.02 sparkGTEx-train.out/parseSalmonReadsSelectCountsCLI.out/onlyNumReads ' > downloadNumReads.sh.batch.02.out 2>&1 &

        setsid sh -c 'time /private/home/aedavids/extraCellularRNA/terra/deseq/bin/downloadNumReads.sh numReadsFileListBatch.txt.03 sparkGTEx-train.out/parseSalmonReadsSelectCountsCLI.out/onlyNumReads ' > downloadNumReads.sh.batch.03.out 2>&1 &
        ```
        
        should easily fit into memory

        ```
        $ ssh mustard
        $ cd /scratch/aedavids
        $ du -sh sparkGTEx-train.out
        212G	sparkGTEx-train.out
        ```

3. use unix paste to create the count matrix

    - create a list of all the part files
    it is very important that the list be sorted. The order must match the
    DESeq colData matrix. given the large number of sample sorting is the only
    way we easily test/debug
    ```
    $ cd /scratch/aedavids
    $ find sparkGTEx-train.out -name "part*.csv" > find.out.listOfNumReadsPartFiles.txt
    $ sort find.out.listOfNumReadsPartFiles.txt > sortedListOfNumReadsPartFiles.txt
    $ head sortedListOfNumReadsPartFiles.txt > testPartFiles.txt
    ```

    run unix paste
    ```
    $ setsid sh -c 'time paste `cat /scratch/aedavids/sortedListOfNumReadsPartFiles.txt` > GTExTrainNumReadsMatrix.tsv' > paste.GTExTrain.out 2>&1 &
    ```

    Sanity test
    ```
    $ ll paste.GTExTrain.out; cat paste.GTExTrain.out
    -rw-r--r-- 1 aedavids prismuser 47 Jan 14 10:13 paste.GTExTrain.out

    real	40m29.662s
    user	27m1.444s
    sys	10m29.626s

    $ wc -l GTExTrainNumReadsMatrix.tsv 
    5387496 GTExTrainNumReadsMatrix.tsv
    ```

    check we have the expecte number of columns. Count the number of tabs in the first line of tsv

    expected number of parts files
    ```
    $ wc -l  /scratch/aedavids/sortedListOfNumReadsPartFiles.txt
    10408 /scratch/aedavids/sortedListOfNumReadsPartFiles.txt
    ```

    we should have 10408 -1 columns
    ```
    $ head -n 1 GTExTrainNumReadsMatrix.tsv | awk -F'\t' '{ print NF-1 }'
    10407

    $ tail -n 1 GTExTrainNumReadsMatrix.tsv | awk -F'\t' '{ print NF-1 }'
    10407
    ```

    add the Name column
    ```
    $ cut -f 1 GTEX-13RTK-1826-SM-5S2P3.quant.sf > transcriptNames.txt
    $ mv GTExTrainNumReadsMatrix.tsv t
    $ paste transcriptNames.txt  t > GTExTrainNumReadsMatrix.tsv

    $ head GTExTrainNumReadsMatrix.tsv | cut -f 1,2,3
    Name	GTEX-1117F-0226-SM-5GZZ7	GTEX-1117F-0526-SM-5EGHJ
    ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|	0.0	0.0
    ENST00000450305.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000002844.2|DDX11L1-201|DDX11L1|632|transcribed_unprocessed_pseudogene|	0.0	0.0
    ENST00000488147.1|ENSG00000227232.5|OTTHUMG00000000958.1|OTTHUMT00000002839.1|WASH7P-201|WASH7P|1351|unprocessed_pseudogene|	764.343	661.807
    ENST00000619216.1|ENSG00000278267.1|-|-|MIR6859-1-201|MIR6859-1|68|miRNA|	0.0	0.0
    ENST00000473358.1|ENSG00000243485.5|OTTHUMG00000000959.2|OTTHUMT00000002840.1|MIR1302-2HG-202|MIR1302-2HG|712|lncRNA|	0.0	0.0
    ENST00000469289.1|ENSG00000243485.5|OTTHUMG00000000959.2|OTTHUMT00000002841.2|MIR1302-2HG-201|MIR1302-2HG|535|lncRNA|	0.0	0.0
    ENST00000607096.1|ENSG00000284332.1|-|-|MIR1302-2-201|MIR1302-2|138|miRNA|	0.0	0.0
    ENST00000417324.1|ENSG00000237613.2|OTTHUMG00000000960.1|OTTHUMT00000002842.1|FAM138A-201|FAM138A|1187|lncRNA|	0.0	0.0
    ENST00000461467.1|ENSG00000237613.2|OTTHUMG00000000960.1|OTTHUMT00000002843.1|FAM138A-202|FAM138A|590|lncRNA|	0.0	0.0
    ```

    copy to gcp project bucket
    ```
    $ gsutil cp GTExTrainNumReadsMatrix.tsv gs://${BUCKET}/quant/sparkGTEx-train.out/
    ```

# failed work arounds
A. transpose the part*.csv files
    use sed to replace new line with comma in our 1 column (it transposes files>
    ```
    setsid sh -c 'time /private/home/aedavids/extraCellularRNA/terra/deseq/bin/createSalmonNumReadsTransposeMatrix.sh sortedListOfNumReadsPartFiles.txt' >> sortedGTExTrainSalmonNumReadsTransposeMatrix.csv 2> createSalmonNumReadsTransposeMatrix.sh.stderr.out
    ```
    
    problem we never found a good way to transpose the file back to the format DESeq expects. I.E. each col is a different sample. We tried using pandas.read_csv().transpose() ran for ever using over 800 GB of memory. eventually terminated. Never completed read. Tried have pandas.read_csv() read in chunks. ran really slowly forced termination after serveral days still had not completed read.

# How the GroupByGeneId count matrices where created
run submitSparkEstimateScalingFactors.sh.  You pass it a txId2GeneId mapping file and a GTEx*NumReadsMatrix.tsv

The GTEx training set was to big to run on spark cluster with 2.8 Tb of memory. We split the train set in to batches by columns and then used paste to put them back together. DO NOT USE SPLIT. it partitions the file by rows. groupBy will not work correctly

next run terra/jupyterNotebooks/cleanUpGeneCountMatrixColNames/ipynb. the groupby column names are of the form  'sum(GTEX-1117F-0226-SM-5GZZ7)' we want 'GTEX-1117F-0226-SM-5GZZ7'

ref: kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv
