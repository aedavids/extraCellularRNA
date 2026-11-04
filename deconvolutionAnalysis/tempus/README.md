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

    ```
    aedavids@mustard $ wc -l mixture.txt 
    137 mixture.txt

    aedavids@mustard $ head -n 1 mixture.txt | tab2newLine | wc -l
    15802

    aedavids@mustard $ head mixture.txt | cut -f 1,2,3
    sampleTitle	GTEX-1117F-0226-SM-5GZZ7	GTEX-1117F-0526-SM-5EGHJ
    A2M	13958.893226398635	57957.76714851816
    ACSL1	3812.7540388156785	756.2916782588303
    ACTA2	8625.168399834383	53611.99355736693
    ACTB	73150.3506185954	79749.48861321727
    ACTG2	1022.7929577581386	6378.947967391779
    ADH1B	66625.29506094292	17860.64439450822
    AHNAK	43081.22905917278	33061.628894290494
    ALB	233.80485221773282	263.02827111982805
    ALDOA	22525.403871789702	33384.09472278025
    ```

