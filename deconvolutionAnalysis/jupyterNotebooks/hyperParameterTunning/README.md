# README.md
Andrew Davidson  
aedavids@ucsc.edu  


date are file system dates. use creation dates in actual file

**Jan 24 09:48 isSignatureMatrixSparce.ipynb**  
  ref: extraCellularRNA/deconvolutionAnalysis/doc/isSignatureMatrixSparce.md
  analysis.bestCuratedGeneConfig.py does not work as well as expected. 
  explores data assumptions for trival case where each category has only 1 gene
  provides a way to explore error and tune algo for larger gene sets
  
***
**Jan 15 10:56 hyperparameterTunningResults4.html**  
**Jan 15 09:29 hyperparameterTunningResults4.ipynb**  

- bake off
- best10CuratedDegree1_ce467ff	mean sensitivity 0.808831  	716	genes


***
**Jan 15 09:29 curatedTunningTargets2.ipynb**  

- best10CuratedDegree1_ce467ff.degree1LUSC10	mean_sensitivity = 0.81 724 genes
- best10CuratedDegree1_ce467ff	                mean_sensitivity = 0.81 716 genes

Jan 15 09:29 curatedTunningTargets.ipynb

Jan 11 10:10 lump.ipynb

***
**Jan 11 10:10 hyperparameterTunningResults3.ipynb**  
**Jan  4 16:07 hyperparameterTunningResults3.2024-01-04.html**  

aedwip 

***
**Jan  4 16:07 hyperparameterTunningResults2.ipynb**  

- rankNGTEx_TCGA

***
**Jan  4 11:26 hyperparameterTunningResults1.ipynb**  
**Dec 28 11:12 hyperparameterTunningResults1.2023-12-28.html**  
**Dec 25 15:37 hyperparameterTunningResults1.2023-12-25.html**  

- includes some explorations. 
- [26] best100Enriched_6_Degree1_selectiveEnrich_LUAD_LUSC_12 best mean sensitivity 0.81 ,  661 genes
- [29] table of sensitivity and specificity of elife types

*** 
**Dec 21 10:51 createDataSets-2023-12-21.html Dec 19 13:38 createDataSets-2023-12-19.html**  

-  best results where mean sensitivity of 0.74 with 969 genes
   metrics for bestNGTEx_TCGA, bestNRemovedGTEx_TCGA, bestNEnrichedGTEx_TCGA, bestNEnrichedDegree1GTEx_TCGA
   * bestNGTEx: select the top N genes for each category
   * bestNRemovedGTEx : genes shared by many categories do not make good biomarkers
   * bestNEnrichedGTEx_TCGA  try to insure each type has a degree1 intersection with a minium of number of genes
   * bestNEnrichedDegree1GTEx_TCGA
    
    
    
** ascending-vs.-DescendingBaseMeanGeneSignatureSelection.ipynb**  
- best10CuratedDegree1_ce467ff has best results. 
- potential bug: sorted base mean in ascending order. this picks potentiall weak signals
- best10CuratedDegree1 is sorted in desecnding order.
- create box plots and histograms of base means to decide if best10CuratedDegree1_ce467ff is valid
