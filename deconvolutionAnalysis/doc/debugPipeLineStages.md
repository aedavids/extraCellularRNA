# Debug PipeLine Stages
Andrew E. Davidson
1/2/24

**pipeline overview**  

ref: deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/hyperparameterTunningResults1.ipynb  

```
best20GTEx_TCGA -> best20RemovedGTEx_TCGA -> best20EnrichedGTEx_TCGA 
    -> best20EnrichedDegree1GTEx_TCGA
                                          
best20GTEx_TCGA -> best20RemovedGTEx_TCGA -> best20Enriched_6_GTEx_TCGA 
    -> best100Enriched_6_Degree1_selectiveEnrich_LUAD_LUSC_{3,6,9,12}
```

** overview of how upstream intersection dictionary is used**

- best20GTEx_TCGA : bestSignatureGeneConfig.py

- best20RemovedGTEx_TCGA : bestRemoveHighDegreeSignatureGeneConfig.py
  + remove genes from upstream intesections with degree > x
  
- best20EnrichedGTEx_TCGA  : EnrichSignatureGeneConfig
  + try to ensure each degree 1 intersection has at least x gene
  +  uses upstream intersection dict to keep track of genes we have already considered. that is stay we are adding "new genes"
  
- best20EnrichedDegree1GTEx_TCGA.sh : byDegreeSignatureGeneConfig.py
  + uses upstream dict to select genes.
  + returns genes in upstream intersections with degree = x
  
- best100Enriched_6_Degree1_selectiveEnrich_LUAD_LUSC_3 : selectiveEnrichSignatureGeneConfig.py
  + ensures each category in categories has at least numberOfGenesToAdd degree 1 genes
  + uses upstream intersection dict to keep track of genes we have already considered. that is stay we are adding "new genes"


**Bug**
file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/hyperparameterTunningResults1.2023-12-28.html

2023-12-28 results where better than 1/2/24 results. The difference was before 1/1/24 we used
best20RemovedGTEx_TCGA.sh used best${topN}GTEx_TCG upstream intersection dictionary. Afte 1/1/24 we swtiched to best${topN}RemovedGTEx_TCGA/training". This intersection dictionary does not have a complete history of genes we have considered. potential we add back poor discriminators. 

/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best20EnrichedGTEx_TCGA/training/best20EnrichedGTEx_TCGA.sh.out/upsetPlot.out

KRT5 was a gene removed upstream, but added back by enrich
```
$ grep 'KRT5' best20_degreeThreshold_10_enrich_3.intersection.dict
    ('ACC', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'BRCA', 'Bladder', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Brain_Substantia_nigra', 'CESC', 'CHOL', 'Cells_Cultured_fibroblasts', 'Colon_Sigmoid', 'Colon_Transverse', 'DLBC', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'GBM', 'HNSC', 'Heart_Atrial_Appendage', 'KIRC', 'KIRP', 'Kidney_Cortex', 'LGG', 'LIHC', 'Liver', 'Lung', 'Muscle_Skeletal', 'Nerve_Tibial', 'OV', 'Ovary', 'PAAD', 'PRAD', 'Pancreas', 'READ', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Thyroid', 'UCEC', 'UCS', 'UVM', 'Uterus', 'Vagina'): [   'KRT5'],
```


