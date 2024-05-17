# add degree 2 genes
Andrew E. Davidson  
aedavids@ucsc.edu  
5/15/24  

ref: 
- [bestCuratedNotes.md](file:///./bestCuratedNotes.md)
- deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/findCandidateEnrichmentBiomarkers.ipynb


**<span style="color:red;background-color:yellow">TODO best10CuratedDegree1_ce467ff</span>**  
what is senstivity of stomic, STAD, ESCA, and GTEx Mucosa? we have some of these values already

```
/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10CuratedDegree1_ce467ff/training/best10CuratedDegree1.sh.out/metrics

$ cut -d , -f 1,6,7 t
id,               specificity, sensitivity
ESCA,             0.999,       0.369
Esophagus_Mucosa, 0.996,       0.991
STAD,             0.999,       0.409
Stomach,          0.999,       0.749
```


## Overview
[hyperparameterTunningResults5.html](file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/hyperparameterTunningResults5.html)
 
over winter break, dec 2023 - jan 2024, I did a lot of deconvolution hyper parmeter tunning. The best results where from run best10LFC_CuratedDegree1. 

specificity is really good  
```
mean_specificity  std_specificity median_specificity numGenes numTypes numDegree1 numAboveThreshold
0.997542	      0.003030	      0.998	             716	  83	   83	      83

LUAD	LUSC	COAD	READ	ESCA	LIHC	STAD	Whole_Blood
1.000	0.994	0.993	0.991	0.998	0.999	1.000	1.000
```
 
elife candidate biomarker sensitivity was not great  
```
mean_sensitivity std_sensitivity median_sensitivity	numGenes numTypes numDegree1 numAboveThreshold
0.786482	     0.206089	     0.833	            716	     83	      83	     62	

LUAD	LUSC	COAD	READ	ESCA	LIHC	STAD	Whole_Blood
0.485	0.595	0.627	0.679	0.396	0.906	0.271	0.987
```


To prepare for new plasma samples I started developing binary Random Forest models on the elife samples.
Our best hyper parameter search results using k-fold validation. "best" was defined as having the highest
AUC, (area under the ROC curve). If we take the best parameter, train a random forest and plot a ROC curve using all the data, the curves look fantastic. Much better the hyper parameter tunning results. I susspect this may because we train using all the data. ie we may be over fitting

I should note we limited the number of biomarkers to 10 to try and minimize overfitting. we only have about 250 elife plasma samples across all classes

```
$ pwd
/private/groups/kimlab/aedavids/elife/hyperparmeterTunning

$ ls *{COAD,Colon,READ}*/*log | tee colorectal.logs.txt
$ ls *{COAD,Colon,READ}*/*csv | tee colorectal.csvs.txt

$ls *csvs.txt
colorectal.csvs.txt  esophagus.csvs.txt  liver.csvs.txt  lung.csvs.txt  stomach.csvs.txt

$ head -n 2 `cat esophagus.csvs.txt ` | cut -d , -f 3,5,7
==> randomForestHyperparmeterEsophagus_ESCASearch.sh.out/randomForestHyperparmeterSearch.csv <==
sensitivity_mean,specificity_mean,auc_mean
0.7861111111111111,0.4238095238095238,0.6049603174603175

==> randomForestHyperparmeterEsophagus_Gastroesophageal_JunctionSearch.sh.out/randomForestHyperparmeterSearch.csv <==
sensitivity_mean,specificity_mean,auc_mean
0.9777777777777779,0.06666666666666667,0.5222222222222221

==> randomForestHyperparmeterEsophagus_MucosaSearch.sh.out/randomForestHyperparmeterSearch.allFloats.csv <==
sensitivity_mean,specificity_mean,auc_mean
0.7861111111111111,0.5142857142857143,0.6501984126984126

==> randomForestHyperparmeterEsophagus_MucosaSearch.sh.out/randomForestHyperparmeterSearch.csv <==
sensitivity_mean,specificity_mean,auc_mean
0.7888888888888889,0.5476190476190477,0.6682539682539683

==> randomForestHyperparmeterEsophagus_MucosaSearch.sh.out/randomForestHyperparmeterSearch.noAUC.csv <==
sensitivity_mean,specificity_mean,max_features
1.0,0.0,1.0

==> randomForestHyperparmeterEsophagus_MuscularisSearch.sh.out/randomForestHyperparmeterSearch.csv <==
sensitivity_mean,specificity_mean,auc_mean
0.6277777777777777,0.5714285714285715,0.5996031746031746
```


our initial london calling nanoporeAdenocarcinomaBinaryClassification.ipynb results using the same biomarkers do not work well. The ROC curves is not good. We are using the nanopore samples as a true hold out. The results are similuar to our k-fold results. a possible explination is that the control samples where from very old patients and not a good match to the age of our ESCA/Adenocarcinoma samples. The 'old' controls where intended for use with alzheimer's plasma samples.

## Try adding biomarkers


file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/intraExtraRNA_POC/jupyterNotebooks/elife/elifeBinaryRandomForestResults.html


The two best elife biomarker sets where  
**GTEx Mucosa**  
```
HUGO_Genes
['FAM135A', 'YOD1', 'THOC3', 'USP6NL', 'NLRX1', 'TXNDC17', 'PGM2', 'DHRS1', 'TTC7B', 'PLD2']

hyperParameterSearchResultsPath:
 /private/groups/kimlab/aedavids/elife/hyperparmeterTunning/randomForestHyperparmeterEsophagus_MucosaSearch.sh.out/randomForestHyperparmeterSearch.csv

 best hyperparameter search results
accuracy_mean accuracy_std sensitivity_mean	sensitivity_std	specificity_mean specificity_std auc_mean  auc_std
0	          0.674286	   0.059446	        0.838889	    0.117063	     0.452381	     0.194015  0.645635	
```

**TCGA ESCA**--
```
HUGO_Genes
['AC012615.3', '(TA)n', 'UBE2SP2', 'HERVFH19-int', 'PRELID1P1', 'LTR106', 'AC010336.3', 'GOLGA8S', 'MER5C', 'CCDC160']

hyperParameterSearchResultsPath:
/private/groups/kimlab/aedavids/elife/hyperparmeterTunning/randomForestHyperparmeterEsophagus_ESCASearch.sh.out/randomForestHyperparmeterSearch.csv

best hyperparameter search results
accuracy_mean accuracy_std	sensitivity_mean sensitivity_std specificity_mean specificity_std auc_mean	auc_std
0	          0.620952	    0.092807	     0.716667	     0.102665	      0.485714	      0.280346	0.60119
```


search decon search resutls for estimate of how to improve these classes

find our estimate of num features notes


### Explore deconvolution upset plots and metrics
ref: deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/findCandidateEnrichmentBiomarkers.ipynb

0) Are there any unused degree genes?

1) Look for degree 2 genes in upset plot we can add.
deconvolutionAnalysis/bin/1vsAll-~gender_category/best10CuratedDegree1.sh
```
/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500LFC_FindAllDegree1_wl500/training/best500LFC_FindAllDegree1_wl500.sh.out/upsetPlot.out/best500LFC_findAllDegree1_wl500.intersection.dict
```



2) Compare Differential Expression Values

3) look at miss classification errors
see metrics

```
/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10CuratedDegree1_ce467ff/training/best10CuratedDegree1.sh.out/metrics

 $ grep 'trueCat\|ESCA\|Esophagus_Mucosa\|STAD\|Stomach' classificationErrors.csv  | cut -d , -f 2,3,4
trueCat,predCat,errorCount
Adipose_Visceral_Omentum,Stomach,1
BRCA,Esophagus_Mucosa,3
CESC,Esophagus_Mucosa,5
COAD,STAD,1
ESCA,BLCA,3
ESCA,Brain_Nucleus_accumbens_basal_ganglia,1
ESCA,Brain_Substantia_nigra,1
ESCA,CESC,2
ESCA,COAD,6
ESCA,Cells_Cultured_fibroblasts,1
ESCA,Colon_Sigmoid,2
ESCA,Colon_Transverse,1
ESCA,Esophagus_Mucosa,11
ESCA,Esophagus_Muscularis,2
ESCA,HNSC,10
ESCA,LUSC,6
ESCA,PAAD,4
ESCA,READ,3
ESCA,STAD,13
ESCA,Skin_Not_Sun_Exposed_Suprapubic,2
ESCA,Stomach,1
ESCA,UCS,1
Esophagus_Mucosa,Adipose_Subcutaneous,1
Esophagus_Mucosa,Esophagus_Gastroesophageal_Junction,1
Esophagus_Mucosa,Minor_Salivary_Gland,1
HNSC,Esophagus_Mucosa,7
KIRP,ESCA,1
LUSC,Esophagus_Mucosa,1
Lung,Stomach,1
Minor_Salivary_Gland,Esophagus_Mucosa,10
OV,ESCA,1
STAD,Adipose_Subcutaneous,2
STAD,Artery_Aorta,1
STAD,Artery_Coronary,3
STAD,BLCA,7
STAD,Brain_Nucleus_accumbens_basal_ganglia,1
STAD,Brain_Substantia_nigra,1
STAD,Breast_Mammary_Tissue,2
STAD,CHOL,2
STAD,COAD,28
STAD,Cells_Cultured_fibroblasts,1
STAD,Colon_Sigmoid,16
STAD,ESCA,17
STAD,Esophagus_Gastroesophageal_Junction,1
STAD,Esophagus_Mucosa,2
STAD,HNSC,2
STAD,LUAD,1
STAD,Lung,2
STAD,PAAD,24
STAD,READ,9
STAD,Stomach,7
STAD,UCS,4
Stomach,Adipose_Subcutaneous,3
Stomach,Adipose_Visceral_Omentum,2
Stomach,Bladder,1
Stomach,Cells_Cultured_fibroblasts,1
Stomach,Colon_Sigmoid,27
Stomach,Esophagus_Gastroesophageal_Junction,19
Stomach,Esophagus_Muscularis,1
Vagina,Esophagus_Mucosa,16
```

adding ESCA and STAD biomarkers may fix 50 miss classification errors.  
adding stomach might fix 70 miss classificaiton errors  
```
trueCat,predCat,errorCount
ESCA,Esophagus_Mucosa,11
ESCA,STAD,13

STAD,ESCA,17
STAD,Esophagus_Mucosa,2
STAD,Stomach,7
```
