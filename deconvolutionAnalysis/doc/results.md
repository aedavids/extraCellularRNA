# Hyperparameter tunning results
aedavids@ucsc.edu 
Dec 20223

## Focus on Lung (tempest is sending us 50 lung cancer samples

### 12/5/23 best20GTEx_TCGA.sh.out/metrics/metricsRounded.csv 
sensitivity is poor  
sensitivity = TP / (TP + FN)  

```
root=/private/groups/kimlab/aedavids/deconvolution/best20GTEx_TCGA/trainingSet/
data=${root}/best20GTEx_TCGA.sh.out/metrics/metricsRounded.csv 
createLungDeconvolutionMetricsCSV.sh $data > data/metrics/best20GTEx_TCGA.Lung.csv 

cut -d , -f 1,5,6,7 best20GTEx_TCGA.Lung.csv
id,         support,specificity,sensitivity
LUAD,       309.0,  1.0,        0.447
LUSC,       301.0,  0.997,      0.495
Lung,       347.0,  0.998,      0.971
Whole_Blood,453.0,  1.0,        0.982
```

**Can we reduce FN?**  

```
cap
cd ~/extraCellularRNA/deconvolutionAnalysis/python
export PYTHONPATH="${PYTHONPATH}:`pwd`"
```

explore LUAD errors  

```
root=/private/groups/kimlab/aedavids/deconvolution/best20GTEx_TCGA/trainingSet/best20GTEx_TCGA.sh.out/
results=${root}/CIBERSORTxFractionsWorkflow.wdl.output/results.txt
expected=${root}/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-20/ciberSortInput/expectedFractions.txt
```

```
cd extraCellularRNA/deconvolutionAnalysis/data/metrics

outDir=exploreClassificationErrors.best20GTEx_TCGA.FN.LUAD.out

python -m analysis.exploreClassificationErrors fn --category LUAD --results $results  --expected $expected --outDir $outDir

python -m analysis.exploreClassificationErrors fp --category LUAD --results $results  --expected $expected --outDir $outDir
```

top LUAD false negatives  
```
outDir=exploreClassificationErrors.best20GTEx_TCGA.FN.LUAD.out

# sort args
# --key use second column
# -n numeric values
# -r revser order
cat falseNegativeGroupByCounts.csv  | sort --field-separator=, --key=2 -n -r | head -n 3
Spleen,77
DLBC,19
Lung,15
```

find genes in the intersection of LUAD and Spleen  
```
python -m analysis.exploreClassificationErrors sg --intersectionDictionary $root/upsetPlot.out/best20.intersection.dict --category LUAD,Spleen --outDir $outDir
```
<span style="color:red">UNEXPECTED RESULTS"</span> There are only three interesections that contain Spleen and LUAD. each of these intersections is composed of many other sets. There are only a hand full of genes shared between LUAD and Spleen. Are the samples mislabled? Quick test: drop these genes?

[   'DES'], [   'KRT6A'], [   'CLU']

60/83 gene signatures have DES in their top 20
47/83 gene signatures have KRTGA in their top 20
25/83 gene signatures have CLU in their top 20

https://www.genecards.org/cgi-bin/carddisp.pl?gene=DES  
https://www.ncbi.nlm.nih.gov/gene/1674#summary  
Summary: This gene encodes a muscle-specific class III intermediate filament. Homopolymers of this protein form a stable intracytoplasmic filamentous network connecting myofibrils to each other and to the plasma membrane. Mutations in this gene are associated with desmin-related myopathy, a familial cardiac and skeletal myopathy (CSM), and with distal myopathies. [provided by RefSeq, Jul 2008]

https://www.genecards.org/cgi-bin/carddisp.pl?gene=KRT6A&keywords=KRT6A  
https://www.ncbi.nlm.nih.gov/gene/3853#summary  
Summary: The protein encoded by this gene is a member of the keratin gene family. The type II cytokeratins consist of basic or neutral proteins which are arranged in pairs of heterotypic keratin chains coexpressed during differentiation of simple and stratified epithelial tissues. As many as six of this type II cytokeratin (KRT6) have been identified; the multiplicity of the genes is attributed to successive gene duplication events. The genes are expressed with family members KRT16 and/or KRT17 in the filiform papillae of the tongue, the stratified epithelial lining of oral mucosa and esophagus, the outer root sheath of hair follicles, and the glandular epithelia. This KRT6 gene in particular encodes the most abundant isoform. Mutations in these genes have been associated with pachyonychia congenita. In addition, peptides from the C-terminal region of the protein have antimicrobial activity against bacterial pathogens. The type II cytokeratins are clustered in a region of chromosome 12q12-q13. [provided by RefSeq, Oct 2014]

https://www.genecards.org/cgi-bin/carddisp.pl?gene=CLU&keywords=CLU  
https://www.ncbi.nlm.nih.gov/gene/1191#summary  
Summary: The protein encoded by this gene is a secreted chaperone that can under some stress conditions also be found in the cell cytosol. It has been suggested to be involved in several basic biological events such as cell death, tumor progression, and neurodegenerative disorders. Alternate splicing results in both coding and non-coding variants.[provided by RefSeq, May 2011]


top LUSC errors  

```
cd extraCellularRNA/deconvolutionAnalysis/data/metrics

root=/private/groups/kimlab/aedavids/deconvolution/best20GTEx_TCGA/trainingSet/best20GTEx_TCGA.sh.out/
results=${root}/CIBERSORTxFractionsWorkflow.wdl.output/results.txt
expected=${root}/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-20/ciberSortInput/expectedFractions.txt
outDir=exploreClassificationErrors.best20GTEx_TCGA.FN.LUSC.out

python -m analysis.exploreClassificationErrors fn --category LUSC --results $results  --expected $expected --outDir $outDir

python -m analysis.exploreClassificationErrors fp --category LUSC --results $results  --expected $expected --outDir $outDir
```


```
cd $outDir
cat falseNegativeGroupByCounts.csv  | sort --field-separator=, --key=2 -n -r | head -n 6
Spleen,48
Lung,12
UCEC,11
HNSC,11
CESC,11
BLCA,9
```

find genes in the intersection of LUSC and Spleen  
```
python -m analysis.exploreClassificationErrors sg --intersectionDictionary $root/upsetPlot.out/best20.intersection.dict --category LUSC,Spleen --outDir $outDir
```
Similar results to LUAD,Spleen. 2 shared genes
[ 'DES', 'CLU']

### 12/6/23 kimlab/aedavids/deconvolution/best20RemovedGTEx_TCGA/training
Removed genes that are poor discriminators. They may be biologically interesting. small improvement in lung sensitivity and specificity. see ~/googleUCSC/kimLab/daniel1:1/2023-12-08/daniel-1-1-2023-12-08.pptx 


**upset plot degree 1 analysis**
several types do not have any uniqu genes. Adding type specific genes should improve deconvolution results

**types that do not have unique genes in our signatue matrix**  
<span style="color:red">There are to many types to pick this off the upset plot. write a tool. the list bellow is incomplete</span>  

Brain_Caudate_basal_ganglia, Brain_Hippocampus, Brain_Amygdala, Brain_Substantia_nigra, Brain_Hypothalamus, Brain_Nucleus_accumbens_basal_ganglia, Brain_Spinal_cord_cervical_c-1, Brain_Putamen_basal_ganglia, Brain_Cortex, Brain_Frontal_Cortex_BA9, Nerve_Tibial

Skin_Sun_Exposed_Lower_leg, Adrenal_Gland, PRAD, KIRP, CESC, Skin_Not_Sun_Exposed_Suprapubic, READ, LGG, Spleen, HNSC, COAD, Esphagus_Muscularis, Adipose_Subcutneeous, Uterus, Liver, SKCM, Adipose_Visceral_Omentum, Colon_sigmoid, Cells_Cultured_fibroblasts, Cells_EBV_transformed_lymphocytes, 

## 12/11/23
using 1vsAll design ~ + geneder + category may out perform ~ + category. 

```
pwd
/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best20GTEx_TCGA/trainingSet/best20GTEx_TCGA.sh.out/metrics

/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/bin/createElifeDeconvolutionMetricsCSV.sh metricsRounded.csv  | cut -d , -f 1,5,6,7
id,support,                                 specificity, sensitivity
COAD,158.0,                                 0.998,       0.437
Colon_Sigmoid,224.0,                        0.987,       0.902
Colon_Transverse,243.0,                     0.997,       0.539
ESCA,111.0,                                 1.0,         0.189
Esophagus_Gastroesophageal_Junction,225.0,  0.996,       0.458
Esophagus_Mucosa,333.0,                     0.998,       0.895
Esophagus_Muscularis,309.0,                 0.995,       0.605
LIHC,223.0,                                 1.0,         0.713
LUAD,309.0,                                 1.0,         0.447
LUSC,301.0,                                 0.997,       0.495
Liver,136.0,                                0.998,       0.934
Lung,347.0,                                 0.998,       0.971
READ,56.0,                                  0.994,       0.821
STAD,225.0,                                 0.998,       0.276
Stomach,215.0,                              1.0,         0.702
```

```
/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category
aedavids@mustard $ cd ../1vsAll-~category/best20GTEx_TCGA/training/best20GTEx_TCGA.sh.out/metrics/

/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/bin/createElifeDeconvolutionMetricsCSV.sh metricsRounded.csv  | cut -d , -f 1,5,6,7
id,support,                                 specificity, sensitivity
COAD,158.0,                                 0.997, 0.272
Colon_Sigmoid,224.0,                        0.986, 0.781
Colon_Transverse,243.0,                     0.995, 0.473
ESCA,111.0,                                 0.996, 0.261
Esophagus_Gastroesophageal_Junction,225.0,  0.995, 0.427
Esophagus_Mucosa,333.0,                     0.993, 0.958
Esophagus_Muscularis,309.0,                 0.995, 0.456
LIHC,223.0,                                 0.998, 0.807
LUAD,309.0,                                 0.998, 0.227
LUSC,301.0,                                 0.992, 0.402
Liver,136.0,                                0.992, 0.779
Lung,347.0,                                 0.998, 0.735
READ,56.0,                                  0.996, 0.5
STAD,225.0,                                 0.995, 0.338
Stomach,215.0,                              1.0,   0.66

```


## 12/14/23

hyper parmater axis:   
best n : increase the number of best genes 
<hr />

hyper parmater axis:   
best n -> remove poor discriminators -> enrich to ensure all types have 3 uique genes  
<hr />
**best20GTEx_TCGA/trainingSet/best20GTEx_TCGA.sh.out/metrics**  
```
createLungDeconvolutionMetricsCSV.sh metricsRounded.csv | cut -d , -f 1,6,7
id,          specificity, sensitivity
LUAD,        1.0,         0.447
LUSC,        0.997,       0.495
Lung,        0.998,       0.971
Whole_Blood, 1.0,         0.982
```

**best20RemovedGTEx_TCGA/training/best20RemovedGTEx_TCGA.sh.out/metrics/ **  
```
createLungDeconvolutionMetricsCSV.sh metricsRounded.csv | cut -d , -f 1,6,7
id,          specificity, sensitivity
LUAD,        0.999,       0.482
LUSC,        0.996,       0.661
Lung,        0.998,       0.977
Whole_Blood, 1.0,         0.958
```


**best20EnrichedGTEx_TCGA/training/best20EnrichedGTEx_TCGA.sh.out/metrics**  
```
createLungDeconvolutionMetricsCSV.sh metricsRounded.csv | cut -d , -f 1,6,7
id,          specificity, sensitivity
LUAD,        1.0,         0.46
LUSC,        0.998,       0.495
Lung,        0.998,       0.974
Whole_Blood, 1.0          ,0.985
```


# 12/15/23 Tunning grid search axis
- change remove degree. we are using 10, try, 5 and 8
- enrich we are using 3. try 1, 5, ...
