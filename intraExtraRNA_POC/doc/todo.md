# Intracellular, Extracellular RNA proof of concept TODO
Andrew E. Davidson
aedavids@ucsc.edu

1. **re-compute "best signatue genes". **
   - select genes that have a base mean > the the average base mean
   - this insures that not only is the gene statitically, it has a strong signal is more likely to be biologically sigifigant 
   

Where was this work done?

a. file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.html
    + This notebook does not do filtering. I think it just aggregates the 1vsAll into a single matrix

b.  make a backup copy of old results
```
# old output
/private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25/ciberSort/signatureGenes.tsv

```

- def CreateBestCibersortGeneSignatureMatrix()
  * I think this might aggregate signatures into a single matrix
  * <span style="color:red">we ran into trouble using long path names with special chars  with cibersort.</span>
  * https://docs.docker.com/storage/volumes/
  * only "[a-zA-Z0-9][a-zA-Z0-9_.-]"
  * best practice use absolute paths
  * consider changing name convention
  * capture baseman threshold

b. /Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.html
    - select "best", "up", "down" gene
    - designed to run on terra "uber workspace" is all the data till on terra?
    - should we create a version that only runs best on mustard?
    - comment out cell [9] it is a stub for 11 def findBestSignatureGenes(deseqDF, signatureGeneConfig):
    - [12] class SignatureGeneConfig
      * args
        + selectGenesOfInterestFunction
        + dataOutputBucketRoot
        + title
        + member funct
          * saveGenesOfInterestToBucketURL
    - tweak [14] def createGTExTCGA_Config_best25()
      * title , may also effect file name
    - [15] <span style="color:red"> this is where the global config obj is set </span>
    - where does data go
      * check bucket "data/1vsAll"
    - [25] moves selected sig genes files to bucket
    

