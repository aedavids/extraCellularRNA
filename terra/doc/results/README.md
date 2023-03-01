# Results
```
Andrew E. Davidson
aedavids@ucsc.edu
Daniel Kim Lab
```

# Table of Contents
- [TODO](todo.md)
- [compute costs](computeCosts.md)
- [clean up](cleanUp.md)
- [../../jupyterNotebooks/cibersort/README.md](../../jupyterNotebooks/cibersort/README.md)

# reference data files
- what salmon ref index did we use?
  * gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.tar.gz
    ```
    ? /private/groups/kimlab/indexes/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx
    $ ls sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/
  complete_ref_lens.bin   info.json         rank.bin             refseq.bin
  ctable.bin              mphf.bin          refAccumLengths.bin  seq.bin
  ctg_offsets.bin         pos.bin           ref_indexing.log     versionInfo.json
  duplicate_clusters.tsv  pre_indexing.log  reflengths.bin
    ```
- /private/groups/kimlab/genomes.annotations/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv


## Cibersort Preliminary results
data in cibersort format was create using .ipynb on mustard. It has been copied over to kimLab/extraCellularRNA/terra/deseq/data to enable upload to https://cibersortx.stanford.edu/upload.php

### deconvolution test 1 best/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:2
Goal: Estimate cibersort docker run time performance. Took about 3.5 days to process GTEX_TCGA training set. (15,801 samples)
   
   
input: extraCellularRNA/terra/deseq/data/private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25/ciberSort:

signatureMatrix has 833 genes: /private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25/ciberSort/signatureGenes.tsv

configuration
    - 100 iterations
    - run mode relative
    
extraCellularRNA/terra/cibersortx/bin/run_cibersortx_fractions.sh.

output files

```
cd /scratch/aedavids/cibersort.out
./ciberSort.debug.out.2022-10-18-07.37.32-PDT 
./2022-09-27-simpleTest/CIBERSORTx_2022-09-27-simpleTest_Results.txt
./GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT/CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt

mv output dir to  /private/groups/kimlab/GTEx_TCGA/cibersort.out/GTEx_TCGA_TrainGroupby_mixture/
```


### results from run Crete (Chania)  sat 10/22/22

```
$ ll /scratch/aedavids/cibersort.out/GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT/CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt 
-rw-r--r-- 1 root root 27M Oct 21 18:42 /scratch/aedavids/cibersort.out/GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT/CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt
```

input was training data set. each sample was from a single type. ie lung, PAAD, ...


approximate total run time 3.5 days
```
  2022-10-18-07.40.54-PDT
-      Oct 21 18:42
-----------------------

 00:20 # 10/18 08:00
+16:00 # 10/19
+24:00 # 10/20
+24:00 # 10/21
+18:42
---------
 83:02

```

### explore results file
15801 samples
```
$ wc -l *Results.txt
15802 CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt
```

number of columns in results
```
$ ls *Results.txt > t
$ sft
$ head -n 1 $f | sed -e 's/\t/\n/g' |wc -l
87
```

These samples are from GTEx. looking at the *Results.txt file alone, it is hard to know what the true tissue type is
Correlation looks good. we expect that samples only deconvole into a single type: TODO write code to check if fraction match expection

```
$ head  $f | cut -f 1,2,3,84,85,86,87
Mixture	ACC	Adipose_Subcutaneous	Whole_Blood	P-value	Correlation	RMSE
GTEX-1117F-0226-SM-5GZZ7	0.00000238818947786	0.32088468153413718	0.00000000000000000	0.00000000000000000	0.98542467644677711	0.92645530582020463
GTEX-1117F-0526-SM-5EGHJ	0.00002175019341662	0.29529977326057638	0.00000000000000000	0.00000000000000000	0.97979080799155893	0.93469518567034882
GTEX-1117F-0726-SM-5GIEN	0.00015259963369839	0.05445776007930397	0.00266438319091741	0.00000000000000000	0.98490561378709340	0.44791559959147931
GTEX-1117F-2826-SM-5GZXL	0.00000000000000000	0.29004114942525649	0.01183538119295482	0.00000000000000000	0.98846407133779690	0.90767508784039519
GTEX-1117F-3226-SM-5N9CT	0.00021702689123660	0.00445898085087070	0.00043860247181520	0.00000000000000000	0.91705227835834124	0.94955257328411946
GTEX-111CU-0326-SM-5GZXO	0.00012839158257415	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.96999619583782004	0.71740096723763391
GTEX-111CU-0426-SM-5GZY1	0.00017911589602714	0.06144797178936864	0.00000000000000000	0.00000000000000000	0.96897428384635731	0.80545057475847193
GTEX-111CU-0526-SM-5EGHK	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.90218261892042595	1.23410866579020806
GTEX-111CU-0626-SM-5EGHL	0.00073656971878856	0.20411159050680064	0.00000000000000000	0.00000000000000000	0.82951499979038523	0.95537827265096409
```

these samples are from TCGA. Looking at the results.txt file alone, it is hard to know what the true cancer type is 
correlation for some of the samples does not look great. Keep in mind we have not done any tunning yet.
```
$ tail  $f | cut -f 1,2,3,84,85,86,87
UVM-WC-A881-TP	0.00000863062071495	0.00309889939958328	0.00000000000000000	0.00000000000000000	0.79120013970914405	0.95817331960717933
UVM-WC-A882-TP	0.00000000000000000	0.00000000000000000	0.01065985495150996	0.00000000000000000	0.64824208271638017	0.95821000145648949
UVM-WC-A883-TP	0.00002268063872719	0.00196036493385387	0.08205200129025236	0.00000000000000000	0.93931661338215500	0.71381242623382180
UVM-WC-A884-TP	0.00015138225222482	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.76022553352495725	0.96080226819233827
UVM-WC-A888-TP	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.82092837581099531	0.96390208810991762
UVM-WC-AA9E-TP	0.00000000000000000	0.00000000000000000	0.00677352970766185	0.00000000000000000	0.70668108537098062	0.95809575383809553
UVM-YZ-A980-TP	0.00095955758606998	0.00314759591266955	0.00000000000000000	0.00000000000000000	0.63927936917453354	0.97024860638613908
UVM-YZ-A983-TP	0.00000894227722533	0.00000000000000000	0.00000000000000000	0.01000000000000001	0.61636688808242612	0.96229386850789200
UVM-YZ-A984-TP	0.00000000000000000	0.00000000000000000	0.00312798052610109	0.00000000000000000	0.79698294355561805	0.95684729590763562
UVM-YZ-A985-TP	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.91638874519350044	0.95127364087000732
```

Our untuned model can only account for 38% of the variation in sample UVM-YZ-A983-TP. UVM is Uveal Melanoma

```
>>> rSquared = 0.616 ** 2
>>> rSquared
0.379456
```

### initial plots

**preliminaryCibersortResults.ipynb**

for each sample/mixture fraction cibersort reports correlaction, p-value, and RMSE

- number of mixtures with correlationThreshold >= 0.9 : 9873 %62.48
- number of mixtures with p-value <= 0.01            : 15154 %95.90
- number of mixtures with correlationThreshold >= 0.9 and pvalue <= 0.01 : 9873
- number of mixtures with correlationThreshold < 0.9 or pvalue > 0.01 : 5928

- whisker plots, suggests cibersort work
  * correlations look high
  * p-values look low
  * not clear what to make of RMSE

I think this is misleading. It does not consider all the information, ie we know what the fractions should be

- Euclidean distance vs. Correlation scatter plot & wisker plot
  * we calculate the distance between the prediced fraction and the expected fraction
  * this is similar to RMSE, it gives us a feel for how close the predictions are to expected
  * whisker plot does not add very much to the story

- tables, column headers colored by data set (GTEx, TCGA)
  * display fraction values across all 83 classes
  * hard to use
  * <span style="color:red">only plotted kidney_cortex (GTEx), traces of kidney cancer.</span>
  * <span style="color:red">Is this noise or indication of early or undiagnosed cancer?</span>
  * <span style="color:red">did a pathologist inspect samples?</span>
  * dropped columns where fraction is very small to make it easier to interpurt table
  * drop columns where col sum is <= 1st quantile ie. 0.01734
['BRCA', 'Brain_Hypothalamus', 'Brain_Substantia_nigra', 'Colon_Transverse', 'ESCA', 'Heart_Atrial_Appendage', 'LGG', 'LIHC', 'LUAD', 'MESO', 'Minor_Salivary_Gland', 'PAAD', 'PCPG', 'PRAD', 'Pituitary', 'SARC', 'Stomach', 'TGCT', 'THCA', 'THYM', 'Testis', 'Vagina']

**contingencyTables.ipynb**
Having 83 different classes makes it hard make sencse of the ciber sort faction by looking at a table. 15,801 rows by 83 columns. One way to make it easier to interpert the data is to use a contigency table. For the Kidney releated classes, we create a indicator variable. For example if a sample as some healthy and KICH signal but not KIRC or KIRB we create a locgical vector [true, true, false, false]. We can then create a table showing how many samples fall into each of the possible permutations. This is a little mis-leading because we do not get any indication about the strength of single components

**volcanoPlots.ipynb**
1vsall plots for kidney releated class and whole blood. TE's are labeled. All plots show lots of TE's

**fractionsAsMulticlassClassification.ipynb**
Cibersort does not use all the avalible information durring train. In our case we know the fractions.
Cibersort fits the fractions to the sample such that RMSE is minimized. We get correlation, RMSE, and p-value. It is hard to determin if fractions are correct. If we treat the fraction as a probablity distribution produced by a multiclassifer, we can call the "class" my simply taking the hard max of the output layer. We can not use standard classifier metrix to evaluate the fractions. **As far as I know this approach has not bee used before**. fractionsAsMulticlassClassification notebook generated a series of confusion matrices. One family of confusion matrices shows count values. The othe family show's row percentages. it for trueth class "A" what percent where TP or type, c, d, e. ...

we have some under reprsented classes so looking at percent is misleading, at same time looking at counts is hard to inperpurted give we have 83 classes.

Preliminary observations:
- 'Adipose_Visceral_Omentum' did not have any TP
- we did not make any predicitons for 'Artery_Tibial' 'Bladder'
- kidney releated class have good accuracy
- some class perform poorly. e.g. artery_aorta, artery_coronary, and artery_tibial


## Software engineering
cibersort took 3.5 days to run on mustard.

1. It is possible to split mixture matrix into parts and run in parralle
   * extraCellularRNA/terra/cibersortx/cibersortParallelization.md
