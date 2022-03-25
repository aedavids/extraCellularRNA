# TODO

## improve spark performance
many of the spark performance problems might be fixed by storing our text/CSV/TSV file in Delta/Parquet format

 Delta (open source file format that uses parquet under the covers). With parquet being columnar many of unused columns can be skipped over vs reading the whole file like csv. 

Here is a link that quickly highlights of parquet.

https://databricks.com/glossary/what-is-parquet

https://medium.com/datalex/5-reasons-to-use-delta-lake-format-on-databricks-d9e76cf3e77d

# To run python unit test from the command line

1. set spark home
```
cd spark-3.1.2-bin-hadoop3.2
export SPARK_HOME=`pwd`
```

2. start extraCellularRNA conda enviroment
```
cd extraCellularRNA
conda activate extraCellularRNA
```

3. modify PYTHONPATH
```
(extraCellularRNA) $ cd extraCellularRNA/terra/deseq/python/test
(extraCellularRNA) $ export PYTHONPATH="${PYTHONPATH}:`pwd`"
```

4. run unit tests
```
(extraCellularRNA) $ cd test
(extraCellularRNA) $ python -m unittest discover .
```

# configure eclipse project
1. in navigator select project properties
2. pyDev - PYTHONPATH
3. select externa libraries tab
4. add source folder
5. enter path to spark-3.1.2-bin-hadoop3.2/python


# required software 
1. [google-cloud-sdk](https://cloud.google.com/sdk/docs/quickstart#installing_the_latest_version)
   - used to move file between Terra workspace and native gcp project and dataproc
2. [apache spark](https://spark.apache.org/)
   - only need if you run unit tests or notebooks locally
   
