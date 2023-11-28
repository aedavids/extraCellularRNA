# TODO
10/4/23

goal build analysis pipeline that will run on phoenix


**status**  

ref: extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
ref: extraCellularRNA/terra/cibersortx/wdl/CIBERSORTxFractionsWorkflow.slurm.sh
ref: https://giwiki.gi.ucsc.edu/index.php/Overview_of_using_Slurm
```
best 2023-Jul-16 21:20:57
best critera was:

    padjThreshold = 0.001
    lfcThreshold = 2.0
    n = 25
    sorted by baseMean     
best 1 vs all results saved to

BUCKET='gs://fc-e15b796f-1abe-4206-ab91-bd58374cc275/data/1vsAll/best-sortedByBaseMean'
${BUCKET}/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25

PRIVATE=/private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/best-sortedByBaseMean
${PRIVATE}/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25
images and intersection.dict

"${BUCKET}/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25GTEx_TCGA_1vsAll-bestN=25-Signature-Genes,-padj-lt-0.001-lfc2-lt=--2.0-or-2.0-lt=-lfc2--min_degree=*.png"

${PRIVATE}/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25GTEx_TCGA_1vsAll-bestN=25-Signature-Genes,-padj-lt-0.001-lfc2-lt=--2.0-or-2.0-lt=-lfc2-.intersection.dict

${PRIVATE}/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25GTEx_TCGA_1vsAll-bestN=25-Signature-Genes,-padj-lt-0.001-lfc2-lt=--2.0-or-2.0-lt=-lfc2--max_degree=*.png
```

**tunning hints**  

there are a couple of class that do not have unique genes. Do they deconvole correctly? Most class have only 1 unique gene
2022-Aug-11 17:18:41 when we sorted by padj instead of base me we can lots of unique genes in each class
there are several genes test that are common to 30 + classes. Does eleminating these genes improve deconvolution results?

## TODO

1. ** create best mixture matrix**
```
CibersortMixtureFactory(
        signatueGeneFilePath,
        groupByGeneCountFilePath, 
        colDataFilePath,
        scalingFactorsPath,
        localCacheDir = 
)
```

```
(extraCellularRNA) aedavids@mustard $ export PYTHONPATH="${PYTHONPATH}:/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python"
(extraCellularRNA) aedavids@mustard $ echo $PYTHONPATH
:/private/home/aedavids/extraCellularRNA/src:/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python
```