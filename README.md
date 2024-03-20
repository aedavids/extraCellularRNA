# extraCellularRNA
Andrew E. Davidson  
aedavids@ucsc.edu  

## find DESeq data sets

anything with de.seq or de-seq is going to be de seq output of either the normalized counts or differential expression variety
 ```
 find /public/groups/kimlab  -name "*de.seq*" -print |grep -v "permission denied"
 ```
 
## conda references
- [https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
- [cheat sheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)


## start the env

```
cd extraCellularRNA
conda activate extraCellularRNA
export PYTHONPATH="${PYTHONPATH}:`pwd`/src"
```

## create the extraCellularRNA environment from yaml file

```
conda env create -f environment.yml
pip install tensorflow
```

## updating dependencies
see [exporting-an-environment-file-across-platforms](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#exporting-an-environment-file-across-platforms)

```
cd extraCellularRNA
conda activate extraCellularRNA
conda env export --from-history > environment.yml
```

## Running Unit test

```
cd extraCellularRNA
conda activate extraCellularRNA
export PYTHONPATH="${PYTHONPATH}:`pwd`/src"
cd src/test
python -m unittest discover .
```

example of how to run a specific test case
```
(extraCellularRNA) $ python plots/test/testUpsetPlots.py  TestUpsetPlots.testFindDegrees
```

## Spark Install
download  spark-3.1.2-bin-hadoop3.2.tgz from https://spark.apache.org/downloads.html

you could probably just install pyspark using pip if you like

## Configuring Visual Studio Code (IDE)
 
[doc/configureVisualStudio.md](./doc/configureVisualCodeStudio.md)


### Warning !!!!!
This branch contains config file for visual studio code. You will need to modify them to match your development enviroment  
    ```
    (extraCellularRNA) $ ls .vscode/
    settings.json  launch.json
    ```

## Deconvolution Hyperparmeter Tunning Overview
**9/8/23 : Status**  
We created some preliminary data that suggests our method may work. This was good enough to advance to candidacy last spring. The analysis was done using a series of notebooks. We will need to clean this up so we can consistently run the analysis using different parmaters, collect the output and meta data in a format that make it easy to track what works and does not work.

**ref**  

1. [advancement proposal](https://docs.google.com/document/d/1I4NWwF2m4UfX-q0JttqMRl8vxfEl-k9h9ziDVbWuRC4/edit?usp=drive_link)
2. [advancement presentation](/Users/andrewdavidson/Documents/UCSC/PhD/advancementTalk/advancementPresentation-final.pptx)

**assumptions**  

1. we ran salmonQuantWorkflow.wdl and 1vsAllTask.wdl on all the GTEx and TCGA sample
2. the results have been aggregated into count matrices (transcripts, and count by genes)
   * data file can be found at /private/groups/kimlab/{GTEx,TCGA,GTEx_TCGA}
   
** preliminary data pipeline **

1. select signatue gene sets
   * extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.html
   * selects set of 1vsAll genes. ie. best n=25, up regulated, down regulated, ...
   * generates upset plots and gene set interections
   * results where copied from GPC bucket to /private/groups/kimlab

2. Create CIBERSORTx mixture matrix
   * extraCellularRNA/terra/jupyterNotebooks/cibersort/createCibersortMixtureMatrix.ipynb
   * scales counts using DESeq estiamted scaling factors
   * output:
     + mixture matrix
     + expected fractions
     + randomized mixture matirx. random shuffle. Does not contain any information we can use this to evaluate how well our model works
     
3. Create CIBERSORTx gene signature matrix
   *  extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
   * given results files extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb and results form 1vsAll analysis, creates a gene signature matrix

4. Run CIBERSORTx
see extraCellularRNA/terra/cibersortx/wdl/README.md for direction on how to run on phoenix/slurm or CLI

** preliminary analysis pipeline**  

extraCellularRNA/terra/jupyterNotebooks/cibersort. See README.md

aedwip test branch merge
