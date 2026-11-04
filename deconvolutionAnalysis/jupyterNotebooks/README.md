# TOC

**explore1vsAllResults.ipynb**  
The basis of our previous biomarker selection algorithms was to select the top N differential express genes with abs(lfc) >=2, padj <= 0.001, sorted by baseMean. This selects biomarkers that have the strongest 'signal'. These biomarkers may not be teh most sepearable feature. This notebook explores an alternative algorithm that selects signifignat genes with the best lfc


**hyperParameterTunning/**  
see readme.md

- debugBest20GTEx_TCGA.sh.out_LUAD_LUSC.ipynb
- fractionsDebug.ipynb
- createSampleTestIntegrationData.ipynb


**signatureMatrixHeatMap.ipynb**
- signature matrix exploration
- creates heat maps of signature matrices

**./randomForestGeneSignatureDeconvolutionPOC.ipynb**  
- trained a random forest model using default parameters on bulk tissue training set
- RepeatedStratifiedKFold() and cross_validate() out performance CIBERSORTx

- <span style="color:red;background-color:yellow">TODO verify results</span>
    * for each class calculate prediction metrics
    * caluclate number of classes with specificity and sensitivity above threshold
    
    


