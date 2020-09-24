

## output file created by ~/aedavids/extraCellularRNA/R/notebooks/kras.ipsc.DESeq.normalize.Rmd
- kras.ipsc.data.bulk_exo.normalized.deseq.biotype.counts.csv

## file created by ~/aedavids/extraCellularRNA/R/notebooks/kras.ipsc_exo_vs_bulk.DESeq.Rmd 
- endings like *{ctrl,kras}*{1,2,3}.biotype.csv 
  * example file  kras.ipsc.data.bulk.data.day.5.ctrl.3.biotype.csv
  * are estimated transcript counts from the quant.sf file produced by salmon.
  * col headers
    +  c("HGNCT", "HGNCG", "BioType", "TPM", "NumReads")


- DESeq.ctrl.sampleType_ex0_vs_bulk.lcfShrink.csv  DESeq.ctrl.sampleType_ex0_vs_bulk.lcfShrink.meta.txt


- DESeq.ctrl.sampleType_ex0_vs_bulk.50.upRegulatedCounts.csv   and down
  * select genes with unshrunk padj < 0.05
  * sort using shrunk lfc
  * most important 50 genes (i.e. most up or down regulated)