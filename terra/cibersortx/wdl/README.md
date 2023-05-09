# CIBERSORTx Fractions WDL workflow overview

```
aedavids@ucsc.edu
3/19/23
```

## Goal:
Our training mixture file is to big to upload to https://cibersortx.stanford.edu/. When we ran the CIBERSORTx fractions docker on mustard it took over 3 days to complete. Our wdl workflow splits the mixture file into parts and uses the scatter/gather pattern to run CIBERSORT on the parts in parrallel 

ref: 
- ../cibersortParallelization.md 
- ../wdlTest/README.md . scatter/gather POC
-  https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather
- ?? we download a zip file from standford that contained a lot of readme files
- you need to register with standford to use the docker. they will send you a token you pass as an argument


getting help:
- https://bioinformatics.stackexchange.com/search?q=wdl
-  https://openwdl.slack.com/ssb/redirect


**Table Of Contents**
- step 1: Create docker image 
- step 2: Create a CIBERSORTx fractions wdl task

## step 1 Create docker image
The cibersortx/fractions from cibersortx.standford.edu can not be called directly from a wdl task. It assume you have mounted your /src/data and /src/datadir to  your local file system. We do not have docker control over mounding. Also the standford docker does not expose a function we can call from our docker task. 

See dockerFile.wdlCibersort

- what is the docker image end point?
  a wdl task is just a bash script. what application does the wdl task need to call?
    ```
    # slack wdl-nightmares 3/14 glennhicky
    
    $ docker pull cibersortx/fractions
Using default tag: latest
latest: Pulling from cibersortx/fractions
Digest: sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc
Status: Image is up to date for cibersortx/fractions:latest
    ```
    
    inspect dumps a lot of information. look for "Entrypoint"
    ```
    aedavids@mustard $ docker inspect cibersortx/fractions
[
    {
        "Id": "sha256:82a17fe6bf93b6a2ae98faefedf04cde66ebc4cdc2e6539a26ac77df353c768a",
        "RepoTags": [
            "cibersortx/fractions:latest"
    ```
    
    our wdl task will need to call
    ```
    "Cmd": [
                "/bin/sh",
                "-c",
                "#(nop) ",
                "ENTRYPOINT [\"./CIBERSORTxFractions\"]"
            ],
    ```
    
    ```
     "Entrypoint": [
                "./CIBERSORTxFractions"
            ],
    ```
***
## step 2: Create a CIBERSORTx fractions wdl task

<span style="color:red">add notes about using docker biuld file in ../wdlTest partition task</span>

validate wdl and run

initialize terminal
```
conda activate extraCellularRNA
export WDL_TOOLS=/private/home/aedavids/extraCellularRNA/java/bin
```

validate
```
java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar validate cibersortxTask.wdl
```

create input json template
```
java -jar /private/home/aedavids/extraCellularRNA/java/bin/womtool-85.jar inputs cibersortxTask.wdl
```


run with test data
```
runCibersortxFractionsTask.sh 2>&1 | tee runCibersortxFractionsTask.sh.out
```

you should find a results file. You should also inspect stdout
```
$ find cromwell-executions -name CIBERSORTx_aedwip_label_Results.txt
cromwell-executions/cibersortxFractionsWorkflow/7a5ee882-dec7-4ef7-95b1-8ab81df79a19/call-cibersortxFractionsTask/execution/CIBERSORTx_aedwip_label_Results.txt
```
results should look something like. cibersort uses a monte-carlo simulation. we do not set the randome seed. Likley your results will be different
```
 cat cromwell-executions/cibersortxFractionsWorkflow/7a5ee882-dec7-4ef7-95b1-8ab81df79a19/call-cibersortxFractionsTask/execution/CIBERSORTx_aedwip_label_Results.txt
Mixture	T1	T2	T3	P-value	Correlation	RMSE
S1	0.99999999999999989	0.00000000000000012	0.00000000000000000	0.00000000000000000	1.00000000000000000	0.17515043817977474
S2	0.00022852600609523	0.99931443200125780	0.00045704199264684	0.00000000000000000	0.99999994770648137	0.17471418357263008
S3	0.00000000000000000	0.00016059614163199	0.99983940385836811	0.00000000000000000	0.99999998280039615	0.33942140553181471
S4	0.50000000110937914	0.49999999889062091	0.00000000000000000	0.00000000000000000	0.97332852678457504	0.34187081755819654
S5	0.49999999967803560	0.00000000060989777	0.49999999971206660	0.00000000000000000	1.00000000000000000	0.51223165638777368
S6	0.00000000000000000	0.49999999971531511	0.50000000028468483	0.00000000000000000	0.96832966378123353	0.35133759608713955
```

```
$ find cromwell-executions -name stdout
cromwell-executions/cibersortxFractionsWorkflow/7a5ee882-dec7-4ef7-95b1-8ab81df79a19/call-cibersortxFractionsTask/execution/stdout
```
