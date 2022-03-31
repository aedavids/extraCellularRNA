# directory overview

The files in this directory are used to test and develop a docker container capable of running DESeq2 from count matrix data.

## 0) Debugging DESeqScript.R
a. configure startRStudioServer.sh to run IMG='aedavids/extra_cellular_rna_2_01'
b. open DESeqScript.R in RStudio
c. uncomment the following lines
```
#
# DEBUG arguments
# allows us to run script in RStudio
#
# getwd()
# countMatrixFile <- "data/1vsAllTest/unitTestGroupByGenesCountMatrix.csv"
# colDataFile <- "data/1vsAllTest/unitTestColData.csv"
# design <- '~ sex + tissue_id'
# referenceLevel <- 'Lung'
# outFile <- 'unitTestDESeqResultOutfile'
# numCores <- 2
# estimateSizeFactorsOutfile <- 'unitTestEstimateSizeFactorsOutfile'
# oneVsAll <- TRUE
# isCSV <- TRUE
```

d. manual execute the lines in DESeqScript.R . <span style="color:red">DO NOT EXECUTE THE CLI PARSE LINES</span>


## 1) setting up the test enviroment on mustard
0. if you need to create a 1vsAll docker image
   follow directions in extraCellularRNA/terra/deseq/bin/dockerFile.1vsAll

1. edit bin/startRServer.sh, configure image and execute
   * IMG='aedavids/extra_cellular_rna_2_01' # starts rstudio-server
   * IMG='aedavids/test-1vs-all-2' # starts rstudio has 1vsall code 
   * IMG='aedavids/test-1vs-all-3-no-rstudio-server' # production candate
2. connect to the container
   ```
   docker exec -it eager_mirzakhani /bin/bash
   ```
3. start another terminal on mustard. use it to edit runner.sh file and exam the R script output files
   ```
   cd /private/home/aedavids/extraCellularRNA/terra/deseq/R
   ```

## 2) run test with mock data
on the docker container
   ```
   # su rstudio
   $ cd /home/rstudio/extraCellularRNA/terra/deseq/R
   $ export PATH=".:${PATH}"
   ```

use runner.sh to driver execute testDESeqScript.R. This hack test script use the mock datafiles.
Note it uses a hack to get the mock data to run. It is not a production quality script.
It is useful for quickly debugging R and DESeq issues

the test data is in extraCellularRNA/terra/deseq/R. The output will be written to runner.sh.out.

Will use mock data files masterCount.tsv  masterColData.tsv 

```
$ runner.sh
# lots of stuff to console. My guess is this is stuff from messages, [R] writes this to stderr

# our print and cat statements wind up in 
$ ls -lt
total 2556
-rw-rw-r-- 1 rstudio rstudio   3391 Jan 25 02:59  DESeqScript.out
```

The data can be found at

```
ls -l masterCount*
masterCountParts.tsv:
total 32M
-rw-r--r-- 1 aedavids prismuser 11M Oct 20 18:41 part-00000-b1374e87-2342-45f3-bb8a-4bfbc588e8ea-c000.csv
-rw-r--r-- 1 aedavids prismuser 11M Oct 20 18:41 part-00001-b1374e87-2342-45f3-bb8a-4bfbc588e8ea-c000.csv
-rw-r--r-- 1 aedavids prismuser 11M Oct 20 18:41 part-00002-b1374e87-2342-45f3-bb8a-4bfbc588e8ea-c000.csv
-rw-r--r-- 1 aedavids prismuser   0 Oct 20 18:41 _SUCCESS

masterCount.tsv:
total 32M
-rw-r--r-- 1 aedavids prismuser 32M Oct 20 18:41 part-00000-d4d920b0-54d0-42bc-bc02-fa9cd1cde034-c000.csv
-rw-r--r-- 1 aedavids prismuser   0 Oct 20 18:41 _SUCCESS
```

you can test the results of a run as follows

```
diff testDESeqScript.expected.out testDESeqScript.out
```



## 2) Running GTEx train 1 vs all tests on mustard
It is faster to debug our R scripts on a mustard rather than using WDL/cromwell or terra. The batch script is the list of command mentioned above. Strange. If you use top you will 1vsAllRunner.batch.sh and its child script 1vsAllRunner.sh
```
$ docker exec --detach --user rstudio friendly_davinci \
    /home/rstudio/extraCellularRNA/terra/deseq/R/1vsAllRunner.batch.sh
```

## 4) testing production script using real data
The production test data Was constructued using following on mustard

```
extraCellularRNA/juypterNotebooks/spark/createTestDESeq2_MasterDataSets.ipynb
extraCellularRNA/terra/deseq/R/createTestColData.Rmd
```

## 5) testing docker container and WDL using cromwell
it will be much faster to debug localling than on terra. See extraCellularRNA/terra/wdl/README.md
