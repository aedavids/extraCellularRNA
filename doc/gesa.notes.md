# Gene Set Enrichment Analysis (GSEA) for interpreting gene expression data

[abstract Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles](https://www.pnas.org/content/102/43/15545.abstract)

ref:

- ~/workSpace/UCSC/BME-230b/hw2/
    * BME-230b-Spring-2019-hw2_question2.ipynb 
        + short cut ~/googleUCSC/kimLab/extraCellularRNA/juypterNotebooks/BME-230b-hw2/BME-230b-Spring-2019-hw2_question2.ipynb
        + for each cluster run gesapy.pre_rank() output /Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/juypterNotebooks/BME-230b-hw2/gseapy.out/prerank_report_GO_Biological_Process_2018_clusterId-0/
        
            - bunch of PNG files and single CSV
            - sample data
                ```
                (base) $ head -n 3 gseapy.prerank.gene_sets.report.csv 
                Term,es,nes,pval,fdr,geneset_size,matched_size,genes

    antigen processing and presentation of peptide antigen via MHC class II (GO:0002495),-0.7851781640420081,-2.098062984733709,0.0,0.0,99,20,"LAG3,FCGR2B,HLA-DOA,HLA-DOB,AP1S2,HLA-DQA2,HLA-DQB1,CTSD,HLA-DQA1,MARCH1,CTSS,IFI30,HLA-DMB,HLA-DRB5,HLA-DPB1,HLA-DPA1,FCER1G,HLA-DRA,HLA-DMA,HLA-DRB1"

    neutrophil mediated immunity (GO:0002446),-0.6469350504674557,-2.3264398631700574,0.0,0.0,488,95,"ANXA2,S100A11,KRT1,CTSA,PGRMC1,S100P,CLEC5A,HSPA1A,SERPINB10,DSP,IL6,PPBP,CXCR2,TXNDC5,CLEC4C,SLCO4C1,HSPA1B,ADGRE5,COTL1,FCGR3B,TNFAIP6,PTX3,CD63,SLC2A3,VCL,CRISPLD2,CLEC4D,FOLR3,TNFRSF1B,ANPEP,MCEMP1,RNASE2,PLAC8,RETN,TLR2,CKAP4,SIRPB1,PYGL,NFAM1,SIRPA,QPCT,NEU1,GSN,CD93,CTSB,OSCAR,HK3,BST1,PTAFR,JUP,TIMP2,LILRB3,ITGAX,C5AR1,FCGR2A,CD33,ITGAM,CTSH,MGST1,CDA,LGALS3,TUBB4B,NPC2,FPR1,LILRB2,PLAUR,HVCN1,SLC11A1,STXBP2,GCA,CLEC12A,PYCARD,CD36,CFP,CTSD,RAB31,S100A12,CD14,GRN,FGL2,CFD,CTSS,PECAM1,CST3,SERPINA1,CD68,FCN1,MNDA,PSAP,CYBB,TYROBP,FCER1G,FGR,ATP6V0C,NME2"
                ```
                
            - rank by the 'nes' column in the csv
                * 'es' col is enrichment score from GESA
                * 'nes' col is normalized enrichment score
                * "The enrichment score is calculated by walking down the list L, increasing a running-sum statistic when we encounter a gene in S and decreasing it when we encounter genes not in S. The magnitude of the increment depends on the correlation of the gene with the phenotype." (Subramanian et al 2005:15546)
    * BME-230b-Spring-2019-hw2_question2.py
- [GENE ONTOLOGY, GO](http://geneontology.org/)
    * world’s largest source of information on the functions of genes. This knowledge is both human-readable and machine-readable
    
## set up
1. Get SCA data into anadata format
2. run PCA, reduce to 50 dimensions
3. nearest neighboors
4. run louvain clustering
5. scanpy plot color by cluster id

## 2.c for each cluster find best matching cell annotation
- this probably works for any given annotation
- uses hypergeometric analysis
  * [hypergeometric distribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution)
  * a discrete probability distribution that describes 
    + the probability of k successes (random draws for which the object drawn has a specified feature) 
    + in n draws, without replacement, 
    + from a finite population of size N 
    + that contains exactly K objects with that feature, 
    + wherein each draw is either a success or a failure. 
  * In contrast, the binomial distribution describes 
    + the probability of k successes in n draws with replacement.
  * multvariate hypergeometric distribution
    + used to calclucate probablity when there are more than 2 classes


Notes from [The hypergeometric test for gene lists] (http://users.unimi.it/marray/2007/material/day4/Lecture7.pdf)

[goStat R package](https://bioconductor.org/packages/release/bioc/html/GOstats.html)
[GOStat vignettes](https://bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsHyperG.pdf)

1. samples are anontated in some way
2. select samples that have a label you are intesrested
3. think of these sample as being white balls and all the rest are black
4. use hypergeometric test

# Data
/public/groups/kimlab/kras.ipsc/bulk.data/day.7/ctrl.1/gencode.salmon.out

```
$ less quant.sf 
Name    Length  EffectiveLength TPM     NumReads
ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|lncRNA|  1657    674.715 0.422396        11.070
ENST00000450305.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000002844.2|DDX11L1-201|DDX11L1|632|transcribed_unprocessed_pseudogene|       632     451.000 0.000000        0.000
```

## Can we use GESA
we have a single bulk sample. clustering genes in a single sample by TPM is probably not useful

can we use hypergometric distribution test to find over or under expressed genes and use these genes to 
form a "gene signature" we can run GESA on to find most likely pathways expressions?

[GSEAPY: Gene Set Enrichment Analysis in Python](https://gseapy.readthedocs.io/en/latest/introduction.html)
- python 
- looks like it it makes it easy to use GESA using pandas in an interactive python console
- I think it only exposes some of the functionality of GESA
    - enrichr, prerank, 
- [gseapy example](https://gseapy.readthedocs.io/en/latest/gseapy_example.html)

[Gene Set Enrichment Analysis (GSEA) User Guide](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html)

[abstract Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles](https://www.pnas.org/content/102/43/15545.abstract)

need to map our gene expression names to Molecular Signature Database, MSigDB. 


## Enrichr
you can upload a list of genes and it will search its data base to find information about the gene, like pathway, ontology, ...

- [enrichr](http://amp.pharm.mssm.edu/Enrichr)
- [you tube tutorial](https://www.youtube.com/watch?v=HfUZdNJ9a3A)


# Crash on court yard

```
2020-07-13 16:58:29,874 Parsing data files for GSEA.............................
2020-07-13 16:58:30,170 Downloading and generating Enrichr library gene sets......
2020-07-13 16:58:30,172 Enrichr library gene sets already downloaded, use local file
2020-07-13 16:58:30,173 User Defined gene sets is given.......continue..........
2020-07-13 18:43:06,383 5103 gene_sets have been filtered out when max_size=500 and min_size=15
2020-07-13 18:43:06,387 No gene sets passed through filtering condition!!!, try new parameters again!
Note: check gene name, gmt file format, or filtering size.
```

```
---------------------------------------------------------------------------
SystemExit                                Traceback (most recent call last)
<timed exec> in <module>

~/miniconda3/envs/extraCellularRNA/lib/python3.7/site-packages/gseapy/gsea.py in prerank(rnk, gene_sets, outdir, pheno_pos, pheno_neg, min_size, max_size, permutation_num, weighted_score_type, ascending, processes, figsize, format, graph_num, no_plot, seed, verbose)
    985                   min_size, max_size, permutation_num, weighted_score_type,
    986                   ascending, processes, figsize, format, graph_num, no_plot, seed, verbose)
--> 987     pre.run()
    988     return pre
    989 

~/miniconda3/envs/extraCellularRNA/lib/python3.7/site-packages/gseapy/gsea.py in run(self)
    458         self._logger.info("Parsing data files for GSEA.............................")
    459         # filtering out gene sets and build gene sets dictionary
--> 460         gmt = self.load_gmt(gene_list=dat2.index.values, gmt=self.gene_sets)
    461 
    462         self._logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))

~/miniconda3/envs/extraCellularRNA/lib/python3.7/site-packages/gseapy/gsea.py in load_gmt(self, gene_list, gmt)
    133             self._logger.error("No gene sets passed through filtering condition!!!, try new parameters again!\n" +\
    134                                "Note: check gene name, gmt file format, or filtering size." )
--> 135             sys.exit(0)
    136 
    137         self._gmtdct=genesets_dict

SystemExit: 0
```
