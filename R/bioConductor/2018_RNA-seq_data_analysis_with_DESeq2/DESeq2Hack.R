# DESeq2 hack
# compare exo cntl 2 exo kras to test our understanding of DESeq2
# TODO switch to R Notebook or Rmd
# ref: 
#   https://github.com/rreggiar/welcome-to-the-kim-lab/wiki/kras.ipsc-data-dictionary
#   https://github.com/rreggiar/welcome-to-the-kim-lab/wiki/Salmon-output-data-dictionary

# tximport loads data into format Bioconductor expects
library("tximport")

#
# 1) find salmon output files, the quant.sf have the transcript ids and counts
# TODO make sure we know what the arguments to Salomon were
# see /public/groups/kimlab/exoRNA-biomarkers-panc/scripts

rootDir='/public/group/kimlab' # root on courtyard and plaz
rootDir='/home/kimlab'         # container root
exoBase='kras.ipsc/exo.data/gen1c.day.5.exo.input.data'
# we only have extra cellular data for day 5
# /public/groups/kimlab/kras.ipsc/exo.data/gen1c.day.5.exo.input.data/{ctrl,kras}.{1,2,3}
# ls /public/groups/kimlab/kras.ipsc/exo.data/gen1c.day.5.exo.input.data/{ctrl,kras}.{1,2,3}/quantFiles/*/quant.sf
# /public/groups/kimlab/kras.ipsc/exo.data/gen1c.day.5.exo.input.data/kras.1/quantFiles/kras.1/quant.sf

createFilePaths <- function(rootDir,  exoBase, label ) {
  # label is either 'kras' or 'ctrl'
  retFiles = c()
  for (i in 1:3) {
    b1 = paste(sep="/", rootDir, exoBase)
    b2 = paste(sep="", label, '.', i, '/quantFiles/', label, '.', i)
    retFiles[i] = paste(sep="/", b1,b2, 'quant.sf')
  }
  
  return( retFiles )
}

exoCntrlFiles = createFilePaths(rootDir,  exoBase, label='ctrl')
exoCntrlFiles
file.exists(exoCntrlFiles)

exoKrasFiles = createFilePaths(rootDir,  exoBase, label='kras')
exoKrasFiles
file.exists(exoKrasFiles)

#
# 2 load the salmon output files
#

# 
# explore the salmon example files using win workshop
# f1 <- "/usr/local/lib/R/site-library/tximportData/extdata/salmon/ERR188356/quant.sf.gz"
# if you do not use txOut you get 'Error in summarizeFail() failed at summarizing to the gene-level'
# workShopSalmon <- tximport(f1, type='salmon', txOut=TRUE)
#
# we get back a list the names "abundance" "counts" 
# "length"  "countsFromAbundance"
# the class of the of the list items is a matrix
# we load 3 salmon files at once then each list item would have 3 matrix. 
# the rows of the matrix are named using the transcript ids
# use rownames() get transcript ids

#
# explore our salmon files
# https://github.com/rreggiar/welcome-to-the-kim-lab/wiki/Salmon-output-data-dictionary
# the name col has multiple gene and transcript ids
# 

txi.exoKrasSalmon = tximport(exoKrasFiles, type="salmon", txIn=TRUE, txOut=TRUE)
# 
# strange. how come no row names on 3 matrix? maybe it some sort of tensor?
#  head( rownames(txi.exoKrasSalmon$counts[,3] ) )
# NULL
#head( rownames(txi.exoKrasSalmon$counts ) )
# [1] "ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|lncRNA|" 

#
# as expected our 'names' encodes multiple names
# we need to create tx2Gene data frame with 2 columns transcriptId, geneId
# figure what what index was used by `$ salmon quant`

# /public/groups/kimlab/kras.ipsc/exo.data/gen1c.day.5.exo.input.data/kras.1/quantFiles/kras.1
# (base) [aedavids@courtyard kras.1]$ cat cmd_info.json |jq
# {
#   "salmon_version": "0.14.1",
#   "index": "/public/groups/kimlab/gen.v31.k.31.v.0.14.index",
#   "libType": "A",
#   "mates1": "output_forward_paired.fq.gz",
#   "mates2": "output_reverse_paired.fq.gz",

# the index file was moved to /public/groups/kimlab/indexes/gen.v31.k.31.v.0.14.index


# salmon idex output
# /public/groups/kimlab/indexes/gen.v31.k.31.v.0.14.index
# (base) [aedavids@courtyard gen.v31.k.31.v.0.14.index]$ cat refInfo.json |jq
# {
#   "ReferenceFiles": [
#     "gencode.v31.transcripts.fa"
#   ]
# }

# 2. create tx2gen file from the reference annotation
# /public/groups/kimlab/genomes.annotations/gencode.31
# has a bunch of csv files. not clear what the source of all these transcript and gene id maps was or how they where constructed
# 
# Roman said "use cut -d '|' " to create the tx2Gene file. "you do not need to touch the quant.sf files" the names are compound how ever also long at the tx2gene values are in the quant.sf files you will be okay

#$ head gen.31.tx.2.gene | cut -d '|' -f 1,2 | sed  's/|/, /'

f <- "/home/kimlab/genomes.annotations/gencode.31/gen.31.tx.2.gene"
gen.31.tx.2.geneDF <- read.table(file=f, header=FALSE, sep ="|")
tx2geneDF <- gen.31.tx.2.geneDF[, c(1,2)]
names(tx2geneDF) <- c("TXNAME", "GENEID")

# tximport
# tximport imports transcript-level estimates from various external 
# software and optionally summarizes abundances, counts, and transcript 
# lengths to the gene-level (default) or outputs transcript-level matrices (see txOut argument).
#
# argument ignoreAfterBar
# whether to split the tx id on the '|' character to facilitate matching 
# with the tx id in tx2gene
#
# txi is the simple list output of the tximport function
# TODO come up with a better name
quantFile = c(exoCntrlFiles, exoKrasFiles)
treatmentFactor <- c( rep('ctrl', length(exoCntrlFiles)), rep('kras', length(exoCntrlFiles)) ) 
treatmentFactor <- factor(treatmentFactor)
txi <- tximport(quantFiles, type="salmon", 
                tx2gene=tx2geneDF, ignoreAfterBar=TRUE)

# > names(txi)
# [1] "abundance"           "counts"              "infReps"             "length"             
# [5] "countsFromAbundance"
# > head(txi$counts)
# [,1]   [,2]   [,3]   [,4]   [,5]   [,6]
# ENSG00000000003.14 35.110 15.711 23.002 35.422 51.369 69.409
# ENSG00000000005.6   1.000  1.000  0.000  0.000  0.000  1.000
# ENSG00000000419.12  5.000  7.000  1.000 17.000  0.000  0.000
# ENSG00000000457.14 53.653 11.778 23.038 18.978 77.632 84.058
# ENSG00000000460.17  2.091  9.149  6.000  5.653  6.178  3.118
# ENSG00000000938.13  2.000  0.000  0.000  0.000  0.000  0.000
# > class(txi$counts)
# [1] "matrix" "array" 
# > 



library("DESeq2")


# 
# samplesMetaDataDF
# must have a least one col
# Rows of correspond to columns of countData
# cols are features that describe the replicants
# the design model use the col as features
samplesMetaDataDF <- data.frame(quantFiles, treatmentFactor)



# TODO come up with better name
# dds is a DESeqDataSet 
# which is a derived classs of RangedSummarizedExperiment 
# which is derived from SummarizedExperiment
dds <- DESeqDataSetFromTximport(txi, colData=samplesMetaDataDF, design = ~ treatmentFactor)

# filter rows: remove genes if they do not have a count of 5 or more in 
# 4 or more samples
keep <- rowSums(counts(dds) >= 5) >= 4
table(keep)
dds <- dds[keep,]

# 
# vst()
# vsd type "DESeqTransform"

# log counts + psudo count of 1 tend to inflate sample variance low counts such
# that it is even larger than biological variation across groups of samples. In
# DESeq2 we therefore provide transformations which produce logscale data such
# that the systematic trends have been removed. Our recommended transformation
# is the variance-stabilizing transformation, or VST, and it can be called with
# the vst

vsd <- vst(dds)

# pca()

# 7.4.3
plotPCA(vsd, "treatmentFactor")
  
# work shop show how to label PCA plot so for a second model feature

#
# 7.5.1 Differtial expression analysis
#

# note this is the vsd
dds <- DESeq(dds)

# res type DESeqResults
res <- results(dds)

# save
write.csv(res, "kras.ipsc.exo.data.gen1c.day.5.exo.de.csv")

# display top genes
head(res[order(res$pvalue),])

# plot counts of top genes
plotCounts(dds, which.min(res$pvalue), "treatmentFactor")

# plot log2 fold change
# makes a so-called "MA-plot", i.e. a scatter plot of log2 fold changes
# (on the y-axis) versus the mean of normalized counts (on the x-axis)
plotMA(res, ylim=c(-5,5))


see page 23 about shrinkage