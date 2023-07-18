# Deconvolution Tunning Notes
```
7/5/2023
Andrew Davidson
aedavids@ucsc.edu
```

background. We originally ran cibersort on mustard in Oct. It took over 3.5 days. There where a couple of mistakes, 

1.  when we selected signature genes we did not set a base mean threshold. This implies we may have selected genes that where statistically signifgant and had a large log fold change how ever had very weak signal. that is to say on a small number of actual transcripts.

2.  When we selected the signature matrix and mixture genes we may not have scaled counts correctly. Turns on we need to divide by the DESeq scaling factors not multiple

I noted these mistakes in my Advancement Proposal on 5/5/23 and my 6/16/23 addendum

## refs:

- /private/groups/kimlab/GTEx_TCGA/cibersort.out/GTEx_TCGA_TrainGroupby_mixture/run_cibersortx_fractions.sh.meta.out  
  + contains docker cli parameters from run that took over 83 hrs.
  + extraCellularRNA/terra/cibersortx/bin/run_cibersortx_fractions.sh
  ```
  $ cat run_cibersortx_fractions.sh.parameters.txt
  docker run --detach --rm -e USERID=30108 
  -v /scratch/aedavids/GTEx_TCGA/geneSignatureProfiles/best/tmp:/src/data 
  -v /scratch/aedavids/cibersort.out/GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT:/src/outdir cibersortx/fractions 
  --username aedavids@ucsc.edu 
  --token 3f561ab6d4cf373d11f23d8e205b4b72 --mixture GTEx_TCGA_TrainGroupby_mixture.txt 
  --sigmatrix signatureGenes.tsv --perm 100 
  --label GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT --QN FALSE --verbose TRUE
  ```
- /private/groups/kimlab/GTEx_TCGA/README.md
  
-  extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
   select sets of genes, and create upset plots
   + run/ran on terra
   + <span style="color:red">7/5/23 about \$1k left, monthly cost about \$153/month port to phoenix</span>
   + findUpRegulatedSignatureGenes() uses baseMean
   + findBestSignatureGenes() cell 9 ???
   + findBestSignatureGenes() cell 11. Does not use baseMean
   + TODO 
     * improve documentation
     * do we want to port this so we can run it on phoenix?
     * make sure we do not over write existing results
     * define findDownRegulated
     * change number change number of selected genes
     * ...
   + results stored in /private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles

- **extraCellularRNA/terra/jupyterNotebooks/cibersort/README.md**
  + cibersort/createCiberSortGeneSignatureMatrix.ipynb
    * bug, we want to fix these one a time to see how they effect results
    * def _createSignatureMatrix(self)
      - divide by scaling factor. do not multiple
      - normalizedDF = transposeGroupByDF *  self.scalingFactors.values

  + cibersort/createCibersortMixtureMatrix.ipynb
    * class CibersortMixtureMatrix
      - _createLabledMixtue(self)
      - bug? divide do not multiple
      - transposeGroupByDF = transposeGroupByDF * self.scalingFactorDF.values
    * createBestMixtureMatrix()
    * createUpMixtureMatrix()
    * TODO create down regulated

## TODO
1. Does dividing by scaling factors improve deconvolution results?
2. Does selecting signature genes using base mean improve deconvolution results?

### 1. Does dividing by scaling factors improve deconvolution results?

**TODO**

- change cibersort/createCiberSortGeneSignatureMatrix.ipynb and cibersort/createCibersortMixtureMatrix.ipynb
  * divide
  * make sure we do not overwrite old buggy values
  * run cibersort using wdl on  phoenix
  * evaluate results
    + extraCellularRNA/terra/jupyterNotebooks/cibersort/
      - preliminaryCibersortResults.ipynb
        * wisker plots
        * fractions table
      - contingencyTables.ipynb
      - volcanoPlots.ipynb
      - fractionsAsMulticlassClassification.ipynb
    + extraCellularRNA/terra/jupyterNotebooks/cibersort/clusters
  
