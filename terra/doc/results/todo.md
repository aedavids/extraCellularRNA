# TODO
```
Andrew E. Davidson
aedavids@ucsc.edu
Daniel Kim Lab
```

## Missing data
we never ran 1vsAll on LAML
- see extraCellularRNA/terra/doc/results/todo-fixLAML.notes.md

## Cibersort run time.
- ~~docker took 3.5 days to run on GTEx_TCGA training set~~
- ~~can we split training sets, run in parallel and test same results~~
- ~~test using 1 100 samples~~
- ~~see extraCellularRNA/terra/cibersortx/README.md for wdl hacks~~

- it may be faster to do hyperpermater tunning on a small subset of our training data set

## GTEX_TCGA 1vsA
- ~~run design ~ gender + category~~
- ~~run design ~ category~~
- ~~download  '~ category' to mustard ~~

## signatureGenesUpsetPlots.ipynb
<span style="color:red">TODO best filter genes with base mean greater than average base mean</span>

- design ~ gender + category
  * select best
  * ~~select up~~
  * select down
  * select smallest possible 'best' (ie try 83 genes)

- design ~  category
  * select best
  * select up
  * select down
  * select smallest possible 'best' (ie try 83 genes)
  
# CiberSort
- should we scale counts (gene signature matrix and mixture matrix)?
  * typically in ml we call data so that count data is in range 0.0 <= count <= 1.0. reduces training time
  * only scale if evidence we get better results
  * UC 1: select gene counts, do not scale
  * UC 2: we have DESeq estiated scaling factors from orginal 1vsAll for training set.
    + we could select genes of interest from non-training data sets and re-calculate DESeq2 estimated scaling facotrs. 
  * UC 3: use min-max scaling
  * UC 4: use estimated scaling factors then min-max
  
## Use DESeq sample scaling factors (normalize for library size and composition)

We use 1vs all  to select best signature genes. When we run cibersort we can to first apply the DESeq scaling factors. for each new data set we first calculate the deseq scaling factors for all the groupby counts. The select the signature genes

## createCiberSortGeneSignatureMatrix
- design ~ geneder + category
  * select best
  * ~~select up~~
  * select down
  * select smallest possible 'best' (ie try 83 genes)

- design ~  category
  * select best
  * select up
  * select down
  * select smallest possible 'best' (ie try 83 genes)
  
- <span style="color:red">tune GEPs. DESeq2 output contains basemean, GEP threshold, baseMean, log2fold, adj-pvalue</span>
  
- create volcano plots of signature genes
  * color ACMG cancer gene sets
    + [amcg 73 genes](https://www.coriell.org/1/NIGMS/Collections/ACMG-73-Genes)

## createCibersortMixtureMatrix.ipynb

- design ~ gender + category
  * select best
  * ~~select up~~
  * select down
  * select smallest possible 'best' (ie try 83 genes)

- design ~  category
  * select best
  * select up
  * select down
  * select smallest possible 'best' (ie try 83 genes)
  
  
- create addtional meta data to select gene sets. ie lncRNA, biotype
/private/groups/kimlab/genomes.annotations/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv

# Ciber sort
- ~~CIBERSORTx-Tutorialipynb~~
  * ~~trival example used to~~
    + ~~figure out how to~~
      - ~~format files~~
      - ~~upload~~
      - ~~interpret plots and results~~
      
- Cibersort GTEX_TCGA design: ~ gender + category
  * best
    + mixture
    + random mixture
        - <span style="color:red">9/15/22 job 10 I think it failed? maybe it is to big to run @ Stanford? </span>
  * up
    + mixture
    + random mixture
 * down
   + mixture
   + random mixture

- tune cibersort
  * running GTEx_TCGA training set took 3.5 days
    - ? can we split training set, run, and merge results?
      + use cromwell, run on mustard
    - does normalizing counts (i.e. [0,1]) improve perfomance?
  * what do these arguments do. do they scale?
    ```
    --rmbatchBmode       <bool>  Run B-mode batch correction [default: FALSE]
    --rmbatchSmode       <bool>  Run S-mode batch correction [default: FALSE]
    --sourceGEPs    <file_name>  Signature matrix GEPs for batch correction [default: 
                          sigmatrix]
    --QN            <bool>  Run quantile normalization [default: FALSE]
    --absolute      <bool>  Run absolute mode [default: FALSE]
    --abs_method    <char>  Pick absolute method ['sig.score' (default) or 
                          'no.sumto1']
    ```

  - ~~explore preliminary results~~ 
    ```
   /scratch/aedavids/cibersort.out
  ./ciberSort.debug.out.2022-10-18-07.37.32-PDT 
  ./2022-09-27-simpleTest/CIBERSORTx_2022-09-27-simpleTest_Results.txt
  ./GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT/CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt  
    ```
- evaluate cibersort using Synthetic

- create addition gene signature meta data
  * indicator one of the 73 AMCG
  * indicator coding non-conding
  * biotypes
  * indicator 
    + /private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.lncRNA.gene.names.txt
      - list of gene names that are in lncRNA
      - has lots of dups
    + TE Annotation meta data public genomes.annotations/genomes.annotations/gencode.39
      - this is meant to merge with TE insertion names (AluYb8, HERVH, etc)
      - ??? ucsc.rmsk.info.txt
      - ??? ucsc.rmsk.insert.info.txt
        + file name means TE insersions
      - columns: TE name, family clade
      - use grep -i Alu

- evaluate cibersort labels
  * ~~for each sample cibersort return a mixture, p-value, and R score. The assumption is that if R is high the fractions are good. This is probably the correct approach when deconvoling unknown samples how ever we are tunning our signature matrix. We know what the signatures, mixtures, and fractions are in advance. Our goal is to tune the signatures to that the fraction match our lables.~~ 
    + can we think about fractions as a discrete probablity distribution. 
      - how can we compare discrete probablity distributions?
      - How can we improve our signature matrix by either adding or removing genes?
        * ?? drop collinar signatures?
        * ~/Documents/books/machinelearningmastery.com/data_preparation_for_machine_learning see chapters 16 and 17 of feature selection and feature importance
    + think of fractions as vectors
      * encode label fractions as 1 hots
      * use ~~confusion matric~~, ROC, and accuracy to evaluate
      * ~~use eclidian distance to evalute~~
      * use cosine similarity
      * for each type calculate residule (i.e. label - fraction)
        - find outliers
        - calucate risdule mean and variance
        - plot histor gram of residual
        - plot box plot of residuals
        
    * cibersort does not use all the information. e.g. we know what the fractions are in advance. How does cibersort perform compare to multi-classifiers like fully connect nerual network, or random forest
- Why are 'nt fractions perferct 1 hots? Could there be traces of cancer?
       
* evaluating off label signales
kidney samples show signs of kidney cancer. Are these signals real or just noise? What method can we use? see contingency tables.

* can we improve model by adding genes for under performing classes or maybe removing some genes?

* more plots
  - ROC
  - louvain clustering and umap visualization
  - bar chart, shows how many samples in each class, identify under reprsented classes
  - re-arange confusion matrix so gtex and releated cancers are close to each other
  - create kidney releated confusion matrix
  
* other method to confirm cibersort factions
  - Josh comp bio SCA home work
    + adj rand index
    + hyper geometric distribution
    + gsea
    + file:///Users/andrewdavidson/workSpace/UCSC/BME-230b/hw2/BME-230b-Spring-2019-hw2_question1.html
    + file:///Users/andrewdavidson/workSpace/UCSC/BME-230b/hw2/euclid_knn.py
    + file:///Users/andrewdavidson/workSpace/UCSC/BME-230b/hw2/BME-230b-Spring-2019-hw2_question2.html
    + file:///Users/andrewdavidson/workSpace/UCSC/BME-230b/hw2/BME-230b-Spring-2019-hw2_question3.html
    + file:///Users/andrewdavidson/workSpace/UCSC/BME-230b/hw2/BME-230b-Spring-2019-hw2_question4.html

# Specific AIM evaluate TE's potential for biomarkers
1. select TE transcripts
2. group by gene
3. run 1vsAll
4. run cibersort

###############AK
click on genom browser
human assempyt hg38
position linc00678 (gene name)
go
scroll down to repeats section , left side, repeatMasker select full and refresh
move mouse down left . in center small text repeating elements by repeatMaster

if we move up to green section LINC00678 shows two isoforms


## are there traces of cancer in health/control GTEx samples?
see kidney table in preliminaryCibersortResults.html. 

```
for each tissue and cancer 
    create table
    sort by tissue,cancer fractions

```
- we expect cancer samples to have some healthy signal
- we did not expect healthy to have cancer signal
  * is this noise? is this early stage diease?
  * disease you die with vs disease you from?
- ? explore patients. older people are at great risk of cancer

### search for metastasis in GTEx subjects

```
for each patient in GTEx
    select all there samples
    plot table
```

## explore TE's
- ~/googleUCSC/kimLab/extraCellularRNA/terra/deseq/python/plots/volcanoPlots.py marks TE red all others black
- te names file was passed on cli googleUCSC/kimLab/extraCellularRNA/terra/deseq/doc/plots/data/gencode.v35.ucsc.rmsk.te.gene.names.txt
```
(base) $ head gencode.v35.ucsc.rmsk.te.gene.names.txt 
L1P5
AluY
L1MB5
AluSc
AluY
HAL1
L2a
L1MA9
L2
L1PA4
(base) $ tail gencode.v35.ucsc.rmsk.te.gene.names.txt
ALR/Alpha
ALR/Alpha
L1PA2
L1PA2
L1PA2
Tigger5b
ALR/Alpha
L1PA3
L1PA3
ALR/Alpha
```


## how evaluate cibersort fractions as if we ran a multiclassifier
- see notebook pages 100, 103

## Can we use Alex's public blood samples
- see notebook pages 103 "1:1 with Daniel"
  * https://elifesciences.org/articles/75181#content "Cancer type classification using plasma cell-free RNAs derived from human and microbes"
  * zotero : kimlab/cancer/alex-datasdet
  * file:///Users/andrewdavidson/Zotero/storage/X6MIFXCU/acc.html
  * plasma cell-free RNAs (cfRNAs) in cancer patients, we profiled cfRNAs in ~300 plasma samples of 5 cancer types (colorectal cancer, stomach cancer, liver cancer, lung cancer, and esophageal cancer) and healthy donors (HDs) with RNA-seq.
- 1. GEO Accession viewer. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174302.
  * https://www.ncbi.nlm.nih.gov/geo/info/overview.html#query
  * search for data https://www.ncbi.nlm.nih.gov/gds/
    GEO DataSets is a study-level database which users can search for studies relevant to their interests. The database stores descriptions of all original submitter-supplied records, as well as curated DataSets. More information about GEO DataSets and how to interpret GEO DataSets results pages can be found on the About GEO DataSets page.
  * https://www.ncbi.nlm.nih.gov/geoprofiles/
    GEO Profiles is a gene-level database which users can search for gene expression profiles relevant to their interests. More information about GEO Profiles and how to interpret GEO Profiles results pages can be found on the About GEO Profiles page.
- Vikas has scripts to download data

# Add our Covid-19 samples

# Data Vis. 

we have 83 types. Hard "see the forest from trees". currated plots ie kidney + kidney cancer do not tell the full story. Potentially other unknow releationships in the data. We want to let the data speak for itself

- see mail folder PhD/DataVisualization (conversation with Chris V)


~~# validate CIBERSORTxFractionsWorkflow.wdl ~~
~~- extraCellularRNA/terra/cibersortx/runCIBERSORTxFractionsWorkflow.sh~~
~~- we ran two tests ~~

~~The resuls should be almost identical. cibersort uses 100 permutation monte-carlo simulation. There does not seem to be a way to set the random seed. Maybe we use numpy.isclose on a couple of rows in the top, middle, and bottom?~~at


# Baseline, time series
"I believe it will be possible to use exosomes to study human biology, disease, and aging. We can test this hypothesis by collecting time series data. We use the sample collected at time zero as our baseline and then look for changes over time. "

Can we get time series data? 
- sample from person with cancer, durring treatment, after then are declared cancer free?
- cell lines? put cells under some sort of stress
- connor is planning to make lung organoids that you can turn on mutant kras in, so in theory you could just track changes over time
- you could try to validate these signal changes using different stage human lung cancer blood samples
