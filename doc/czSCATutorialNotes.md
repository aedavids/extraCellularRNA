
# 1. Data Processing
Seurat is an R package designed for QC, analysis, and exploration of single cell RNA-seq data.
Bioconductor is a open-source, open-development R project for the analysis of high-throughput genomics data, including packages for the analysis of single-cell data.
Scanpy is a Python package similar to Seurat

##  AnnData Object
- AnnData.X stores the count matrix
- AnnData.obs stores metadata about the observations (cells)
- AnnData.var stores metadata about the variables (genes)
- AnnData.uns stores any additional, unstructured information we decide to attach later

to construct an AnnData object from data frames use 
```
adata = sc.AnnData(X = count_dataframe, obs = metadata_dataframe)
```

to save to disk
```
adata.write('../data/brain_raw.h5ad') ## the h5ad extension is AnnData-specific
```

An [RNA spike-in](https://en.wikipedia.org/wiki/RNA_spike-in#:~:text=An%20RNA%20spike%2Din%20is,known%20as%20a%20control%20probe.)
is an RNA transcript of known sequence and quantity used to calibrate measurements in RNA hybridization assays, 
such as DNA microarray experiments, RT-qPCR, and RNA-Seq. A spike-in is designed to bind to a DNA molecule 
with a matching sequence, known as a control probe

 
# 2. quality control

Since there is currently no standard method for performing scRNAseq, the expected values for the various QC
measures that will be presented here can vary substantially from experiment to experiment. Thus, to perform
QC we will be looking for cells which are outliers with respect to the rest of the dataset rather than 
comparing to independent quality standards. Consequently, care should be taken when comparing quality metrics
across datasets collected using different protocols.
 
 
 see scanpy pp.calculate_qc_metrics(). It returns two dataframes: one containing quality control metrics about cells,
 and one containing metrics about genes.
 
```
 qc = sc.pp.calculate_qc_metrics(adata, qc_vars = ['ERCC'])# this returns a tuple of (cell_qc_dataframe, gene_qc_dataframe)
                                 # ask for the percentage of reads from spike ins
                                
cell_qc_dataframe = qc[0]
gene_qc_dataframe = qc[1]
```

### control for library size
remove genes with low read counts


### detect genes
make sure that the reads are distributed across the transcriptome. Thus, we count the total number of unique 
genes detected in each sample. i.e. use a histogram. remove cells, samples, with low gene counts


### spike ins
Another measure of cell quality is the ratio between ERCC spike-in RNAs and endogenous RNAs. 
This ratio can be used to estimate the total amount of RNA in the captured cells. Cells with a high 
level of spike-in RNAs had low starting amounts of RNA, likely due to the cell being dead or stressed 
which may result in the RNA being degraded.


### filter
use scanpy pp.filter_cells() 

### quailty control
It is typically a good idea to remove genes whose expression level is considered "undetectable". 
We define a gene as detectable if at least two cells contain more than 5 reads from the gene. 
However, the threshold strongly depends on the sequencing depth. It is important to keep in mind that 
genes must be filtered after cell filtering since some genes may only be detected in poor quality cells.

### normalize count

# analysis

### differential expression
[differential expression](https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/05-diffexp.html)

* there are several approaches to this task:
    - Look for upregulation of marker genes for cell types of interest (compared to the rest of the dataset)
    - Compare the complete gene expression profiles between groups
    - Use automated methods to compare cells of interest to databases of cell type expression profiles to 
      combine clustering and annotation. what is an example of this?
        
For well-defined cell types, we expect marker genes to show large differences in expression between the cell type of 
interest and the rest of the dataset, allowing us to use simple methods. 

Important note! For differential expression, we need to use the raw values stored in adata.raw. 
With differential expression, we want to account for both the center and spread of the expression in each group. 
Recall that when we normalized our values, we standardized the distribution of each gene across cells to be 
centered at 0 and scaled with variance 1. So, when calculating differential expression, we should use the 
raw values (post-QC, pre-normalization).
