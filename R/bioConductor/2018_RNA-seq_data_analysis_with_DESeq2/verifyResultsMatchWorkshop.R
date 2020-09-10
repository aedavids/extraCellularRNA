# make sure we know how data was prepared for section 7.5
# there was a lot of prepocessing in early sections

# load the DESeq version of a sumarized experiment data set
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# do not retain genes if they do not have a count of 5 or more
# in 4 or more samples. these genes do not have no statistical power 
# to detect differences, and no information to compute distances 
# between samples.
keep <- rowSums(counts(dds) >= 5) >= 4
dds <- dds[keep,]

# normal way to run DESeq2
dds <- DESeq(dds)
# console output
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
res <- results(dds)


# display top genes
head(res[order(res$pvalue),])
