# Deconvolution Analysis Pipeline Overview

Andrew E. Davidson  
aedavids@ucsc.edu  

goal: user has selected a set of genes of interest and wants to run cibersort and stats on a particular data set

Select of genes set of interest is a manual process

**ref**
- extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
    * Dose not have good documentation
    * runs on terra
    * loads 1vsAll results (DESeq2 results file)
        + there is a file for each GTEx/TCGA tissue type
    * defines functions that select sets of genes that are of interest. E.G upregulated, ...
    * for each 1vsAll Result file
        + selects genes of interest DESeq results and saves them to a file
    * creates 
        + upset plots and intersection.dict

- extraCellularRNA/terra/jupyterNotebooks/cibersort/scalingBug/createCibersortMixtureMatrix.ipynb
    * good documentation. (overview, input, output)

**input**
- experiment id
- path to gene set of interest
- path to train/validate/test dataset
- path to 1vsAll results
- path to 1vsAll estimated scaling factors
- natural languge description
    * notes that will be include in meta data

**algo**
- create signature matrix
- create mixture Matrix
    * ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/scalingBug/createCibersortMixtureMatrix.ipynb
    * input:
        1. a list of signature genes
        2. a gene count matrix
        * The row ids are gene name, the column names are the sample ids
        3. a DESeq ColData matrix.
        * contains sample meta data
        4. DESeq estimated scaling factors
        * adjust each sample to account for libaray size and library composition
    * run cibersort
    * run stats and error metrics
    * save/store
    * create meta data file

# 10/5 design notes
I think we will get better resuse by creating bigger modules. It is hard to  parameterize slurm scripts.We want consistency. ie. run end to end when ever we get a new set of candidate biomarkers. We want to make it easy to re-run, capture pedigree, output, ...

should get better reuse. run from cli, easy to use with slurm, easier to create docker and WDL

1. create a python class to run/produce cibersort input file, upset plots, interesection dictionaries
2. for each hypothesis create a python driver script
    * no CLI. hard code all the arguments
    * input counts, meta data, ...
    * define a filter class
3. create slurm script
    * hard code the driver script