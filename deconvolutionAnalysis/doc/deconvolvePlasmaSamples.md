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

** <span style="color:red;background-color:yellow">best10CuratedDegree1_ce467ff  signatureGenes.tsv. We do not need to calcualte a new signature matrix</span>**
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

- upstreamPipeline.py 
  * /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/pipeline/upstreamPipeline.py

- CibersortMixtureFactory.py
  * /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/cibersortMixtureMatrixFactory.py
    + <span style="color:red;background-color:yellow">has a main()</span>  
    +  <span style="color:red;background-color:yellow">How do we scale the plasma samples?</span>
       - elife counts where already normalized. 
       - for new esoph samples complete seq will normalize counts.
       - typically in machine learn you learn your feature scaling factors on your training set. You apply the same factors on your hold out. You do not recalculate scaling factors on hold out set. It is assume the hold out is drawn from the same population as the training set. Our plasma samples are not from the same population as our training samples. Our hypothesis is that the biomarkers we discover using bulk tissue will work in plasma how ever the distributions , ie the weights, will be different. <span style="color:red;background-color:yellow">lets just see what happens. ie we will not recacluation normalization </span>. Our assumption is the mixture will still be a linear combination of the bulk signature matrix. **PlanB** can we use low rank matrix factorization to calculate the weight of the signature matrix?
       ![Marked logo](file:////Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/deconvolutionAnalysis/doc/img/plasmaLRMF_signatureMatrix.png)
         * <span style="color:red;background-color:yellow">create dumby scaling factor. [1,1,1,...] so we can reuse code</span>
   + colData is used to create 'labels'
   + <span style="color:red;background-color:yellow"> need to map gencode 37 -> 39</span> see intraExtraRNA_POC/python/src/models/randomForestHyperparmeterSearch.py, and elifeUtilities.py
   + _select()
     - will select the normalized subset of genes that comprise the signature matrix
   + _createLabeledMixtue()

**step 3: run deconvolution**  
aedwip

**step 4: analysis**  
aedwip

