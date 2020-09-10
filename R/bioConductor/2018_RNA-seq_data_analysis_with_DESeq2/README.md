# RNA-seq workflow: gene-level exploratory analysis and differential expression

see zotero PDF

ref:
- [workshop 2018 RNA-seq data analysis with DESeq2](https://bioconductor.github.io/BiocWorkshops/rna-seq-data-analysis-with-deseq2.html)
- [bioconductor 2018 workshops](https://bioconductor.github.io/BiocWorkshops/)
- [this might be old? http://master.bioconductor.org/packages/release/workflows/html/rnaseqGene.html](http://master.bioconductor.org/packages/release/workflows/html/rnaseqGene.html)
- [docker containers for bioconductor not sure if they will work with 2018 workshop](https://hub.docker.com/r/bioconductor/bioconductor_docker)
- [biomanager, do not use normal R install()](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
- [install](http://master.bioconductor.org/install/#install-bioconductor-packages)


# install 
use:
- bioconductor/bioconductor_docker:devel
- download source but does not compile or install
```
BiocManager::install("Bioconductor/BiocWorkshops")
```

install the packages we need for tutorial. need use refresh for pkss to  show up in the packages tab in RStudio
```
use example : BiocManager::install(c("DESeq2", "tximport"))

DESeq2
tximport
apeglm
AnnotationHub
ReportingTools
Glimma
splatter
zinbwave
readr

check for
library("DESeq2")
library("tximport")
library("AnnotationHub")
library("ReportingTools")
library("Glimma")
library("splatter")
library(zinbwave)

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("jsonlite")
library("readr")
library("airway")
library("ggplot2")

```


## saved image aedavids/biocworkshop2018desq2

