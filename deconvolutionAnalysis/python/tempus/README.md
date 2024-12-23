Andrew E. Davidson  
<aedavids@ucsc.edu>  
12/22/24  

Copyright (c) 2020-2023, Regents of the University of California All rights reserved. <https://polyformproject.org/licenses/noncommercial/1.0.0>  

code to run deconvolution on tempus dilution samples. We want to run different gene signature matrices.

# Design overview

**reverse eng the deconvolution hyperparameter tunning pipeline**  

best10GTEx_TCGA.sh
  
/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/bin/pipeline.sh  
python -m pipeline.upstreamPipeline

/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/pipeline/upstreamPipeline.py  

* _step1()
  * /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/driver.py
    * runSelectGenesOfInterest()
      * for each deseq result filePath
        * retDict[fileName] = signatureGeneConfig.findGenes(deseqDF, fileName)

**<span style="color:red">TODO: tempus psudo code use runSelectGenesOfInterest()</span>**  
    * The list will only have a single results file annotated_norm_counts.csv

/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/analysis/bestSignatureGeneConfig.py  

    * select using padj, lfc, basemean
    * return top n

upstreamPipeline.py  _step3() create signature matrix  
/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/cibersortSignatureMatrixFactory.py  
/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/test/testCibersortSignatureMatrixFactory.py  

**annotated_norm_counts.csv format**  

    ```bash
    $ head annotated_norm_counts.csv  | cut -d , -f 1,2,3,4 | sed -e "s/,/, /g"
    gene_id, gene_name, gene_biotype, SLDK3_T1_100_S1_L007
    (A)n, (A)n, Microsatellite, 0.5450089567650692
    (AAA)n, (AAA)n, Microsatellite, 0
    (AAAAAAC)n, (AAAAAAC)n, Microsatellite, 0
    (AAAAAAG)n, (AAAAAAG)n, Microsatellite, 0
    (AAAAAAT)n, (AAAAAAT)n, Microsatellite, 0
    (AAAAAC)n, (AAAAAC)n, Microsatellite, 0
    (AAAAACA)n, (AAAAACA)n, Microsatellite, 0
    (AAAAACC)n, (AAAAACC)n, Microsatellite, 0
    (AAAAACT)n, (AAAAACT)n, Microsatellite, 0
    ```

* self._createsignatureGeneList()
* _read_ColData()
* read groupby gene count data ie. annotated_norm_counts.csv
  * remove the  gene_name, gene_biotype columns

**example of sample gene signature format**  
/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10GTEx_TCGA/training/best10GTEx_TCGA.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-10/ciberSortInput  

    ```bash
    aedavids@mustard $ head -n 1 signatureGenes.tsv | tab2newLine | wc -l
    84

    aedavids@mustard $ head -n 1 signatureGenes.tsv | tab2newLine | head -3
    name
    ACC
    Adipose_Subcutaneous

    aedavids@mustard $ wc -l signatureGenes.tsv 
    137 signatureGenes.tsv
    aedavids@mustard $ head signatureGenes.tsv | cut -f 1,2,3,4
    name ACC Adipose_Subcutaneous Adipose_Visceral_Omentum
    A2M 62731.79337360358 56168.65579716901 51871.17804397286
    ACSL1 4999.261152193253 34138.35718602691 35610.426049294525
    ACTA2 6693.980752374152 36679.40152547475 20359.410963253456
    ACTB 205730.14963986518 145117.6246826885 123237.7319031974
    ACTG2 187.3809586420608 5772.336256960208 2250.2555131980876
    ADH1B 4442.995421528372 117999.9270660066 75394.94568758561
    AHNAK 14044.409164631541 125952.97911245028 74832.62840703489
    ALB 410.4155640992933 228.55648904836545 489.0402528089585
    ALDOA 116311.7735104082 26562.637044600327 22489.62560500435
    ```

* _createSignatureMatrix(self)
  * transposeGroupByDF = self.groupByCountDF.transpose(copy=True)

    ```python
    # join the colData, we need the 'category' col so we can
    # calculate the  signature gene mean value for each category 
    joinDF =  pd.merge(left=normalizedDF, 
                        right=self.colDataDF.loc[:,["sample_id", "category"]], 
                        how='inner', 
                        left_index=True, 
                        right_on="sample_id")      
    self.logger.debug(f'joinDF:\n{joinDF}')

    # calculate the expected values for each category      
    genesDF = joinDF.loc[ :,self.geneListsorted + ["category"] ] 
    if self.useMedian :
        # weird duplicated log so I can set debugger break points
        self.logger.info(f'useMedian : {self.useMedian} calling median()')
        signatureDF = genesDF.groupby("category").median()
    else:
        self.logger.info(f'useMedian : {self.useMedian} calling mean()')
        signatureDF = genesDF.groupby("category").mean()
    
    # convert to cibersort expected upload format
    ciberSortSignatueDF = signatureDF.transpose(copy=True)
    ciberSortSignatueDF.index.name = "name"
    
    self.ciberSortSignatueDF = ciberSortSignatueDF
    # weird. cciberSortSignatueDF.columns.name = 'category'. This name is not
    # saved by pd.to_csv(). Set to none to make it easier to write unit test
    self.ciberSortSignatueDF.columns.name = None
    ```

