# directory overview

The files in this directory are used to test and develop a docker container capable of running DESeq2 from count matrix data.

## 1) setting up the test enviroment on mustard
1. bin/startRServer.sh
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
   $ export PATH=".:${PATH}
   ```

in rstudio you will need to install the argparse project

use runner.sh to driver execute testDESeqScript.R. This hack test script use the mock datafiles.
Note it uses a hack to get the mock data to run. It is not a production quality script.
It is useful for quickly debugging R and DESeq issues

you can test the results of a run as follows

```
diff testDESeqScript.expected.out testDESeqScript.out
```

## 3) testing production script using reall data
The production test data Was constructued using following on mustard

```
extraCellularRNA/juypterNotebooks/spark/createTestDESeq2_MasterDataSets.ipynb
extraCellularRNA/terra/deseq/R/createTestColData.Rmd
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



