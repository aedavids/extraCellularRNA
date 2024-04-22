# How to deconvolve plasma samples
Andew E. Davidson  
aedavids@ucsc.edu  
4/22/24  


Reverse engineer deconvolution hyperparameter tunning pipeline.


## Overview

best hyperparameter results: /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/bin/1vsAll-~gender_category/best10CuratedDegree1_ce467ff.sh  

**step 1: constuct a signature matrix from best10CuratedDegree1_ce467ff results**

- 1-vs-all results: 
  * /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10CuratedDegree1_ce467ff/training/best10CuratedDegree1.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-10

- upstreamPipeline.py 
  * /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/pipeline/upstreamPipeline.py
  
- cibersortSignatureMatrixFactory.py
  * deconvolutionAnalysis/python/pipeline/dataFactory/cibersortSignatureMatrixFactory.py
    + <span style="color:red;background-color:yellow">has a main()</span>
    + _createSignatureMatrix()
      - scales counts
      - colData is used to calculate the signature gene mean value for each category 

** <span style="color:red;background-color:yellow">best10CuratedDegree1_ce467ff  signatureGenes.tsv</span>**
``` 
/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10CuratedDegree1_ce467ff/training/best10CuratedDegree1.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-10/ciberSortInput/signatureGenes.tsv

$ wc -l signatureGenes.tsv 
717 signatureGenes.tsv

head -n1 signatureGenes.tsv | tab2newLine | wc -l
84

 $ head -n1 signatureGenes.tsv | cut -f 1,2,3,83,84
name	ACC	Adipose_Subcutaneous	Vagina	Whole_Blood
```

**step 2: construct mixture matrix from plasma samples**  
aedwip

**step 3: run deconvolution**  
aedwip

**step 4: analysis**  
aedwip

