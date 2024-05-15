# add degree 2 genes
Andrew E. Davidson  
aedavids@ucsc.edu  
5/15/24  

## Overview
 hyperparameterTunningResults5.html](file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/hyperparameterTunningResults5.html)
 
over winter break, dec 2023 - jan 2024, I did a lot of deconvolution hyper parmeter tunning. The best results where from run best10LFC_CuratedDegree1. 

specificity is really good  
```
mean_specificity  std_specificity median_specificity numGenes numTypes numDegree1 numAboveThreshold
0.997542	      0.003030	      0.998	             716	  83	   83	      83
```
 
sensitivity was not great  
```
mean_sensitivity std_sensitivity median_sensitivity	numGenes numTypes numDegree1 numAboveThreshold
0.786482	     0.206089	     0.833	            716	     83	      83	     62	
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
