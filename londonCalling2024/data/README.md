# Meta data prep notes

create a tidy meta data file for /private/groups/kimlab/vikas/nanopore/promethion/barretts/analysis/complete-seq/R_results /private/groups/kimlab/aedavids/londonCalling2024/data/

** 1. save meta data to csv files**   
the 'batch' worksheet have the meta data for each of the nanopore runs

for each batch work sheet in spread sheet I 
    - save as -> csv
    - you will be prompted to save only selected worksheet
```
$ ls *.csv
barretts_sample_details_batch4.csv  barretts_sample_details_batch2.csv
barretts_sample_details_batch3.csv  barretts_sample_details_batch1.csv
```

**2. paste a column on each csv file identifying the batch**  
we can use this col to debug

```
$ cat addBatchCol.sh
#!/bin/bash

set +x
batchNumber=$1

# count the number of lines
# we use '<' so that wc does not print file name
l=`wc -l < *batch${batchNumber}.csv`


colOut=batchColFile.txt
'rm' ${colOut}
for i in `seq 0 $l`;
do
    echo "batch${batchNumber}" >> ${colOut}
done

paste $colOut *batch${batchNumber}.csv
```

```
(base) $ addBatchCol.sh 1 > colData1.csv
(base) $ addBatchCol.sh 2 > colData2.csv
(base) $ addBatchCol.sh 3 > colData3.csv
$ addBatchCol.sh 4 > colData4.csv
```

**3. concatenate all file into 1 files **  
``` 
$ cat  colData1.csv  colData2.csv  colData3.csv  colData4.csv > colDataIncomplete.csv
```

use an editor to remove duplicate header lines
```
batch2	sample,sex,age,race,histology,barcode
batch3	sample,sex,age,race,histology,barcode
batch4	sample,sex,age,race,histology,barcode
```

use an editor to change header from
```
 $ head -n 1 colDataIncomplete.csv
batch1,﻿sample,sex,age,race,histology,barcod
```

to
```
 $ head -n 1 colDataIncomplete.csv
batchId,﻿sample,sex,age,race,histology,barcod
```

**4. convert file format to unix**  

```
dos2unix colDataIncomplete.csv > t
mv t colDataIncomplete.csv
```


```
scp * mustard:/private/groups/kimlab/aedavids/londonCalling2024/data
```

**5. use jupyter notebook  to add sampleId col that matches normalized_counts_for_andy.csv**  
extraCellularRNA/londonCalling2024/jupyterNotebooks/createTidyMetaData.ipynb

