# README.md

**Table of Contents**
<span style="color:red">incomplete</span>

notebooks are listed in chronologic order 

terra/jupyterNotebooks/
- createDESeqData.ipynb 11/17/21
  * GTEx
  * creates Train/validate/test data sets and colData

- createDESeqDataTCGA-terra.ipynb 
  * prototype 
  * create colData train/validate/test data sets
  * does not save results

- createGroupByGenesCounts.ipynb 

- TCGA-createGroupByGeneDESeqDataSets

- a) generateTCGAMatrixCreationScripts.ipynb 4/27/22
  * run on mac
  * creates script used to transfer salmon quant files to jupyter notebook vm running in terra
  * finds results from 'salmonTarQuantWorkflow v 4'
    + TCGA workspace data varies from cohort/study group to studyGroup
    + bam files may be in tar, or gzip, contain paired or single end reads, and include replicants
    + It was tricky and expensive to sort this out. for each data model there maybe multiple 
    columns contains salmon quant files. This notebook identifies the correct files
 * sample output for TCGA_READ
   ```
   $ pwd
   /Users/andrewdavidson/googleUCSC/kimLab/terraDataModels
   (base) $ ls test-aedavids-proj/TCGA/TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab/generateTCGAMatrixCreationScripts.ipynb.out
   TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_quantFile.csv
   TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_105_copyFromTerraToLocal.sh
   TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_105_copyFromTerraToNativeGCP.sh*
   TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_colData.csv
   
   $ head -n 2 TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_105_copyFromTerraToLocal.sh
   mkdir -p ./TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab/
   gsutil -m cp gs://fc-secure-8a69fc00-b6c9-4179-aee5-f1e47a4475dd/34b2bbfb-4f9a-41d4-bfd8-b55a8e1987de/quantify/ef52a514-2fc0-4e85-946a-2bbbbc56ab96/call-salmon_paired_reads/READ-AF-2687-TP.quant.sf.gz ./TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab/
   
   ```
  
- b) createCountMatrix.ipynb 
  * terra uber workspace file date 5/27/22
  * uses generated localization scripts created by generateTCGAMatrixCreationScripts.ipynb to copy salmon files to local VM
  * use unix cut and past to create tsv file from the numReads column
  * tcga output gs://fc-e15b796f-1abe-4206-ab91-bd58374cc275/data/matrices/NumReads/
  
- generateTerraCleanUpScripts.ipynb
  creates scripts remove deprecated salmon quant and aux files from the TCGA data set
  

## chronological ordering to jupyter notebook transformations
Terra does a good job processing large batches of samples in parallel. 
When you want to aggregates teh output of your pipeline into a single data structure.
It is typically easier to use juypter notebooks that runs on a single large machine.
A problem with jupyter notebook is that the transformations are typically done in 
several steps, each in a separate notebook. It is easy to loose track of the order
of execution. 

To recover the order of exection I used the git commit. The ordering is not perferct
because work was done on a branch then merge -squash back to main

```
$ git ls-files | while read file; do git log -n 1 --pretty="$file, %ad" -- $file; done  | tee t
```

sort t in excel by date

```
cleanUpGTExDataSets.ipynb,                          Fri Mar 25 13:31:37 2022 -0700
cleanUpGeneCountMatrixColNames.ipynb,               Fri Mar 25 13:31:37 2022 -0700
cleanUpTerraDataSets.ipynb,                         Fri Mar 25 13:31:37 2022 -0700
countQuantFiles.ipynb,                              Fri Mar 25 13:31:37 2022 -0700
createDESeqData.ipynb,                              Fri Mar 25 13:31:37 2022 -0700
createFailedSampleDataSet.ipynb,                    Fri Mar 25 13:31:37 2022 -0700
createGTExTissueDataSets.ipynb,                     Fri Mar 25 13:31:37 2022 -0700
createParseSalmonReadsSplitBatchs.ipynb,            Fri Mar 25 13:31:37 2022 -0700
howToUse_gsutil.ipynb,                              Fri Mar 25 13:31:37 2022 -0700
salmonQuantWorkflowCostAnalysis.ipynb,              Fri Mar 25 13:31:37 2022 -0700
sparkBlockMatrixTransform.ipynb,                    Fri Mar 25 13:31:37 2022 -0700
splitIntoTissueSampleSets.ipynb,                    Fri Mar 25 13:31:37 2022 -0700
uploadTerraDataModelTest.ipynb,                     Fri Mar 25 13:31:37 2022 -0700

cibersort/CreateGeneSignatureMatrixOverview.ipynb,  Fri Sep 9 15:35:13 2022 -0700
cibersort/createCiberSortGeneSignatureMatrix.ipynb, Fri Sep 9 18:31:44 2022 -0700

cibersort/mixture.txt,                              Sun Sep 11 10:11:02 2022 -0700
cibersort/signatureGenes.txt,                       Sun Sep 11 10:11:02 2022 -0700

TCGA-createDESeqDataSets.ipynb,               Thu Sep 1 12:24:31 2022 -0700
TCGA-createGroupByGeneDESeqDataSets.ipynb,    Thu Sep 1 12:24:31 2022 -0700
cibersort/CIBERSORTx-Tutorial.ipynb,          Thu Sep 1 12:24:31 2022 -0700
createCountMatrix.ipynb,                      Thu Sep 1 12:24:31 2022 -0700
createDESeqDataTCGA-terra.ipynb,              Thu Sep 1 12:24:31 2022 -0700
createFailedSampleDataSet-TCGA.ipynb,         Thu Sep 1 12:24:31 2022 -0700
createGTEx-TCGADESeqDataSets.ipynb,           Thu Sep 1 12:24:31 2022 -0700
exploreGTExDESeqColData.ipynb,                Thu Sep 1 12:24:31 2022 -0700
exploreUpsetPlotInteresections.ipynb,         Thu Sep 1 12:24:31 2022 -0700
generateTCGAMatrixCreationScripts.ipynb,      Thu Sep 1 12:24:31 2022 -0700
signatureGenesUpsetPlots.ipynb,               Thu Sep 1 12:24:31 2022 -0700

cibersort/createCibersortMixtureMatrix.ipynb, Tue Sep 6 15:55:28 2022 -0700
```

* signatureGenesUpsetPlots.ipynb
  - originaly written to run on terra
  - selects set of 1vsAll genes. ie. best n=25, up regulated, down regulated, ...
  - generates upset plots and gene set interections
  - results where copied from GCP bucket to /private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles
  
   * cibersort/createCibersortMixtureMatrix.ipynb
   * scales counts using DESeq estiamted scaling factors
   * output:
     + mixture matrix
     + expected fractions
     + randomized mixture matirx. random shuffle. Does not contain any information we can use this to evaluate how well our model works
  
 
