
# bulk deconvolution sensivity is not great
Andrew E. Davidson   
aedavids@ucsc.edu  
3/31/24  

Daniel is meeting with colaborators next week. Do we have preliminary evidence our method works for ESCA?

**Can we just use Esophagus_Mucosa bio markers? I.e. health vs cancer?**  

extraCellularRNA/deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/hyperparameterTunningResults5.html

![hyperparameterTunningResults5 sensitivity and specificity table](file:./esophSensitivitySpecificitityTable.png)

**did any runs perform well for ESCA?***
Yes. maybe we can manually add genes from best1CuratedDegree1_6457b56 results.
```
find . -name metricsRounded.csv | tee t

$ head -n 1 /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10GTEx_TCGA/training/best10GTEx_TCGA.sh.out/metrics/metricsRounded.csv 
id,precision,recall,f1-score,support,specificity,sensitivity,tp,fn,fp,tn

$ grep ESCA $f | cut -d , -f 7 | sort -d | tail -n 4
0.631
0.685
0.685
0.775


$ grep ESCA $f | cut -d , -f 1,7 | grep '0.775\|0.685\|0.685\|0.631'
./best1CuratedDegree1_6457b56/training/best1CuratedDegree1.sh.out/metrics/metricsRounded.csv:ESCA,0.775
./best5CuratedDegree1_ce467ff/training/best5CuratedDegree1.sh.out/metrics/metricsRounded.csv:ESCA,0.685
./best3CuratedDegree1_de50c43/training/best3CuratedDegree1.sh.out/metrics/metricsRounded.csv:ESCA,0.685
./best2CuratedDegree1_6457b56/training/best2CuratedDegree1.sh.out/metrics/metricsRounded.csv:ESCA,0.631
```

# do our intracellular biomarkers work on elife ?
ref: extraCellularRNA/intraExtraRNA_POC/jupyterNotebooks/elife/lungCancer/randomForestIntraCellularLungCancerBiomarkersOnExtracellularSamples.ipynb

1. load data
2. how many esoph and healthy control samples do we have?

## <span style="color:red;background-color:yellow">TODO</span>
- consider using ESCA genes from best1CuratedDegree1_6457b56
- intraExtraRNA.elifeUtilities.loadElifeLungTrainingData fix features argument
  * ~~rename we made it generic, maybe wrapper it it new func is generic, Lung call new generic~~
  * check yNP. should health control = zero?

