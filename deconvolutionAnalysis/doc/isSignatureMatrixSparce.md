# Is Signature matrix sparce? Is there a bug?
Andrew E. Davidson  
aedavids@ucsc.edu  
1/15/24  

In general we get mean sensitivity of 0.81. When we added a few extra LUSC degree 1 genes mean mean sensitivity did not change much how ever LUSC sensitivity dropped

**?Bug 1?**
best10Curated has 716 genes, we expect 830

we select degree 1 genes that are unique in window lenth 1000 of best genes. ie lfc, padj, sorted by base mean in decending order

**What should we expect the signature matrix to looks like?** 

- **<span style="color:green">check a category with high specificity</span>**

- **<span style="color:green">verify we find degree 1 genes correctly</span>**
use upsetPlot to create intersection dictionary for top 1000 genes. Do not train. Search for ABHDF in the intersection dictionary

- Look at ABHDF values in training data to verify 1vsAll results. 

- Could there be a sorting bug in the code that creates signature and mixture matrix? when we saved the estimated scaling factors we did not save the gene ids. 

**We made an assumption. Is it valid?**  

- we select degree 1 genes found in a window of length 1000 best genes ranked by baseMean
- it is likely many categories will express these genes at some level
- lets say one of our bio markers is "g1" and lfc >= 2.0
- we assume the baseMean of g1 > the expressed value in all other categories
  + <span style="color:red">FALSE</span>
    * see 'What is rank of of ABHD5 in Uterus?" bellow
    * check of Best500GTEx_TCGA find ABHD5 in only in Uterus_vs_all.results

- TODO: AEDWIP
It is possible that another category express g1 with a higher base mean outside of the 1000 gene window. ie there where 1000 other genes that had a higer base mean. This seem unlikely given that deseq normalization account for library size and composition

- TODO : AEDWIP
Are we calculating the signature matrix correclty? 
    + I do not think we should use base mean
    + I think we need to do something like
      ```
      sigDF['luad', 'g1'] = countDF['luad', g1].mean()
      sigDF['liver', 'g1'] = countDF['liver', g1].mean()
      sigDF['Whole_Blood', 'g1'] = countDF['Whole_Blood', g1].mean()
      ```
   + I think we can estimate the gene signature
   ```
   lfc = x/baseMean; geneSigMatrx[i,j] = x 
   ```

## quick check explore best1CuratedDegree1

best1CuratedDegree1 is simple. There are only 83 genes in the signature matrix. 

```
cd /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best1CuratedDegree1/training/best1CuratedDegree1.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-1/ciberSortInput
```

Genes names are the row index. They are in alphabetic order. 
```
 $ wc -l signatureGenes.tsv 
84 signatureGenes.tsv

# explore row idx
aedavids@mustard $ cut -f 1 signatureGenes.tsv | head -n 3
name
ABHD5
AC022149.1

aedavids@mustard $ cut -f 1 signatureGenes.tsv | tail -n 3
ZBTB4
ZNF395
ZNF516

# explore 1st header row
$ head -n 1 signatureGenes.tsv | cut -f 1,2,3
name	ACC	Adipose_Subcutaneous

aedavids@mustard $ head -n 1 signatureGenes.tsv | cut -f 82,83,84
Uterus	Vagina	Whole_Blood
```

**Which category was ABHD5 a degree 1 bio marker?**  

ABHDF is degree 1 genes for Uterus. It has a negative lfc. The signature matrix value is 453.19
```
cd /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best1CuratedDegree1/training/best1CuratedDegree1.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-1

$ grep ABHD5 *.results
Uterus_vs_all.results:
name,  baseMean,         log2FoldChange,    lfcSE,             stat,             pvalue,                padj
ABHD5, 1424.84191261374, -2.07498350972574, 0.119833002337334, -17.315626490644, 3.5859817507590704e-67, 3.23471258948033e-65

# Uterus is col 82.

$ head -n 2 signatureGenes.tsv | cut -f 1,81,82,83,84
name	UVM	                Uterus	          Vagina	         Whole_Blood
ABHD5	417.37879461437257  453.1932664458746 1188.512451206417 1494.6709524356854
```

**What is the Rank of ABHD5 in Uterus?**  
13 out 83. on the small size  
```
# take the first row of data 'ABHD5'
head -n 2 signatureGenes.tsv | tail -n 1 >t

# edit t. remove the gene name

# convert tabs to new line
tab2newLine t > tt

# numeric sort
aedavids@mustard $ sort -n tt | > ttt

$ grep -n ^45 ttt
13:453.1932664458746
```

**What categories have smaller values for ABHD5 in the signature matrix?**  

```
# deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/isSignatureMatrixSparce.ipynb
categories with signature matrix less the ABHD5, Utereus degree 1 gene
Brain_Amygdala                           335.922066
Brain_Anterior_cingulate_cortex_BA24     386.884412
Brain_Caudate_basal_ganglia              386.563951
Brain_Cortex                             362.353517
Brain_Frontal_Cortex_BA9                 380.226687
Brain_Hippocampus                        328.161104
Brain_Nucleus_accumbens_basal_ganglia    406.307070
Brain_Putamen_basal_ganglia              339.124982
Liver                                    365.776347
Ovary                                    430.322131
UVM                                      417.378795
```

** explore deseq results for ABHD5**  

<span style="background-color:yellow">deseq results are sorted by padj</span>  
The estimnated scaling factors are in the same row order as original count matrix

```
cd /private/groups/kimlab/GTEx_TCGA/1vsAll
grep -n ABHD5  *results > grep.ABHD5.results.out

# python see isSignatureMatrixSparce.ipynb
df.sort_values(by="log2FoldChange", ascending=True)

	category	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
24	Cells_EBV*	1424.841913	-2.195980	0.107318	-20.462444	4.654597e-93	3.471387e-92
78	Uterus	    1424.841913	-2.074984	0.119833	-17.315626	3.585982e-67	3.234713e-65
54	Ovary	    1424.841913	-2.014958	0.106562	-18.908773	9.657828e-80	4.181927e-78
79	UVM	        1424.841913	-1.845651	0.158850	-11.618825	3.306522e-31	3.460344e-30

1	Adipose*	1424.841913	0.852824	0.061331	13.905320	5.880345e-44	3.678126e-43
65	Skin_not	1424.841913	0.966015	0.057970	16.664160	2.388198e-62	1.362160e-61
81	Whole_Blood	1424.841913	2.995298	0.042534	70.420976	0.000000e+00	0.000000e+00
```

**<span style="background-color:yellow;color:red">How do we reconsile sig matrix with deseq results</span>**  

- Why is ABHD5 a degree 1 genes associated with Uterus? 
  * Cells_EBV-transformed_lymphocytes is lower
    + TODO.  test
    + we only took 1 degree 1 genes.
    + did  Cells_EBV and Uterus have other degree1 genes with a higher baseMean? Even if they did by definition they would be out outside of our 1000 gene window
    + did these other classes really have 1000 "better signature genes?"
  
  * Whole has highest lcf value
  * <span style="background-color:yellow;color:red">maybe we need to rank using a combintation of baseMean and lfc?</span>
  try to get average count lfc = count / base mean. sort by this. what happens if 
      + bm high and lfc high 
      + bm high and lfc low (negative)
      + bm low and lfc high
      + bm low and lfc low (negative)
      


**<span style="color:red">What can we say about degree one genes?</span>**  
- ABHD5 is a degree 1 genes for Uterus
- it's lfc is negative how ever it is not the smallest value the signature matrix.
- the brain class are problematic

# estimate signature matrix values
Is this formula correct?

```
Uterus_vs_all.results: ABHD5

baseMean = 1424.84191261374
lfc =  -2.07498350972574 

lfc = count / baseMean
avg count = 2** -2.07498350972574 * 1424.84191261374 = 338.16947615327723
```











<hr />
## Bone Yard
<hr />

**find the category with the biggest absolute value of ABHD5**  
```
# take the first row of data 'ABHD5'
head -n 2 signatureGenes.tsv | tail -n 1 >t

# edit t. remove the gene name

# convert tabs to new line
tab2newLine t > tt

# numeric sort
aedavids@mustard $ sort -n tt |head

328.16110387461237
335.9220664266857
9819857462
362.3535174009732
365.77634663886533
380.2266870682288
386.56395139379237
386.88441168980563
406.3070704017742

aedavids@mustard $ sort -n tt |tail
2888.293579363596
3086.5277426171997
3136.072152104618
3425.0553831459706
3569.8169248709614
3944.28945042935
4292.063167002956
4831.956592432293
5181.799428125866
8337.443626048012
```

** which category was ABHD5 a degree 1 gene?**  

```
cd /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best1CuratedDegree1/training/best1CuratedDegree1.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-1

$ grep ABHD5 *.results
Uterus_vs_all.results:
name,  baseMean,         log2FoldChange,    lfcSE,             stat,             pvalue,                padj
ABHD5, 1424.84191261374, -2.07498350972574, 0.119833002337334, -17.315626490644, 3.5859817507590704e-67, 3.23471258948033e-65

# Uterus is the last col 82. col=geneName
$ head -n 1 signatureGenes.tsv | cut -f 81,82,83,84
UVM	Uterus	Vagina	Whole_Blood

$ head -n 2 signatureGenes.tsv | tail -n 1  | cut -f 81,82,83
417.37879461437257	453.1932664458746	1188.512451206417
```

**<span style="color:red">ABHDF is clearly not the best biomarker for Uterus</span>**  
It looks like ABHD5 is a good bio marker for SARC?
```
# 8337 is col 66
$ head -n 2 signatureGenes.tsv | tail -n 1  | tab2newLine |grep -n 8337
66:8337.443626048012

 $ head -n 2 signatureGenes.tsv  | cut -f 1,66
name	SARC
ABHD5	8337.443626048012
```

DESeq results for ABHD5 in SARC
```
/private/groups/kimlab/GTEx_TCGA/1vsAll

aedavids@razzmatazz $ grep -n  ABHD5 SARC_vs_all.results
34790:
name,    baseMean,         log2FoldChange,      lfcSE,              stat,              pvalue,             padj
"ABHD5", 1424.84191261374, -0.0629790662352483, 0.0892637414737385, -0.705539171845903, 0.480474763527371, 0.567170744664102

```

# estimate signature matrix values
Is this formula correct?

339.125 is at index 4 in cell sorted list of ABHD5 values in gene signature matrix
```
Uterus_vs_all.results: ABHD5

baseMean = 1424.84191261374
lfc =  -2.07498350972574 

lfc = count / baseMean
count =  2**-2.075 * 1424.842 = 338.166
avg count = 2** -2.07498350972574 * 1424.84191261374 = 338.16947615327723
```

1363.545 is at index 51 in cell sorted list of ABHD5 values in gene signature matrix
```
SARC_vs_all.results: ABHD5
baseMean =  1424.84191261374
lfc =  -0.0629790662352483

lfc = count / baseMean
average count = 2** -0.0629790662352483 * 1424.84191261374 = 1363.980286841443
```
