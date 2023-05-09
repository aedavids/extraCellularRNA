
#### AEDWIP 05/09/23
It would be best to re-run 1vs all if we want to publish the entire 77K differential expressed genes. In our actual work we only use maybe the best 25. We can just run 1 vs all on the LAML set and filter

# missing LAML in signature matrix


**summary**
When we lumped all the TCGA data sets together we missed LAML. If we want the LAML samples in our 1vsAll we would need to reconstuct teh GTEX_TCGA groupby data sets, upload them to terra and re-run the 1vsAll wdl 

## figure out where we/how we create the signature matrix

**file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/terra/jupyterNotebooks/cibersort/CreateGeneSignatureMatrixOverview.html**

- POC shows how to create data sets for CIBERSORT

**file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.html**

- create cibersort data files

**file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.html**

- looks like a notebook that runs on terra
- findUpRegulatedSignatureGenes()
- findBestSignatureGenes() empty function!
- findBestSignatureGenes


# uber workspace, did we run 1vs all LAML?
https://app.terra.bio/#workspaces/test-aedavids-proj/uber/data

LAML is not present in either the uber workspace tables  'GTEx_TCGA_1vsAll' or  'design_b_category_GTEx_TCGA_1vsAll'

There are two table in the uber terra workspace

1. 'GTEx_TCGA_1vsAll'  https://app.terra.bio/#workspaces/test-aedavids-proj/uber/data

   * these results where downloaded to /private/groups/kimlab/GTEx_TCGA/1vsAll
   * design : ~gender + category
   ```
   aedavids@mustard $ gsutil ls gs://fc-e15b796f-1abe-4206-ab91-bd58374cc275/f568dca6-4a0a-4de6-ae70-dcd4c189ae5e/deseq_one_vs_all/b194ea58-32a0-4507-8eae-f130541c580c/call-one_vs_all/Adipose_Subcutaneous_vs_all.results
   
   aedavids@mustard $ head Adipose_Subcutaneous_vs_all.results 
mean of normalized counts for all samples 
log2 fold change (MLE): category Adipose Subcutaneous vs not Adipose Subcutaneous 
standard error: category Adipose Subcutaneous vs not Adipose Subcutaneous 
Wald statistic: category Adipose Subcutaneous vs not Adipose Subcutaneous 
Wald test p-value: category Adipose Subcutaneous vs not Adipose Subcutaneous 
BH adjusted p-values 
design: ~gender + category 
"name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
"(GGT)n",63.499961398417,-5.35234066006956,0.140879048789305,-37.9924531437912,0,0
"AC004816.1",69.4869446421764,-5.78142850669159,0.121901507332856,-47.4270469101355,0,0
   ```
2. 'design_b_category_GTEx_TCGA_1vsAll' https://app.terra.bio/#workspaces/test-aedavids-proj/uber/data

   * results downloaded to /private/groups/kimlab/GTEx_TCGA/1vsAll-~category
   * design: ~category
   ```
   aedavids@mustard $ gsutil cp gs://fc-e15b796f-1abe-4206-ab91-bd58374cc275/submissions/92f62668-c42c-4efb-b642-43308d24d93c/deseq_one_vs_all/d7f49088-7cd9-4436-90eb-20ecf62c5cc1/call-one_vs_all/Pituitary_vs_all.results .
   
   aedavids@mustard $ head Pituitary_vs_all.results 
mean of normalized counts for all samples 
log2 fold change (MLE): category Pituitary vs not Pituitary 
standard error: category Pituitary vs not Pituitary 
Wald statistic: category Pituitary vs not Pituitary 
Wald test p-value: category Pituitary vs not Pituitary 
BH adjusted p-values 
design: ~category 
"name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
"(AAGG)n",2.3825527148895,4.33213966335474,0.114430100234554,37.8583926298667,0,0
"AC010618.1",0.55757499953234,2.91187383123218,0.0553357498231346,52.6219277869945,0,0
   ```
   
# LAML samples are  not in the  TCGA data set

There are 33 cohort types in coldata.csv
```
aedavids@mustard $ pwd
/private/groups/kimlab/TCGA/colData
aedavids@mustard $ ls TCGA*colData.csv | wc -l
33
```

How ever the the training data sets combined colData.csv files only has 32 cohorts
```
aedavids@mustard $ pwd
/private/groups/kimlab/TCGA/colData/trainingDataSets
aedavids@mustard $  cut -d , -f 4 TCGA-TestColData.csv  | sort | uniq | grep -v Cohort | wc -l
32
```

LAML samples ids are not in /private/groups/kimlab/TCGA/trainingDataSets

```
head -n 1 GTEx_TCGA_*Groupby.csv | comma2newLine > $d/sampleNames.txt

$ for sid in `cat lamlSampleIdx.txt`;
> do
> echo $sid
> grep $sid sampleNames.txt 
> done
```




