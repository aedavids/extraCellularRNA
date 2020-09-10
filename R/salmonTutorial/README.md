# Salmon tutorial
ref [https://combine-lab.github.io/salmon/getting_started/](https://combine-lab.github.io/salmon/getting_started/)

shell scripts dl_tut_reads.sh*  quant_tut_samples.sh:

sample data is plant

reads are paired. i.e. left/right 1/2

"Once you have your quantification results you can use them for downstream analysis with differential expression tools like swish, DESeq2, edgeR, limma, or sleuth. Using the tximport package, you can import salmon’s transcript-level quantifications and optionally aggregate them to the gene level for gene-level differential expression analysis. You can read more about how to import salmon’s results into DESeq2 by reading the tximport section of the excellent DESeq2 vignette. "


"When you are building a salmon index, please do not build the index on the genome of the organism whose transcripts you want to quantify, this is almost certainly not want you want to do and will not provide you with meaningful results."
### AEDWIP
build the index using reads from the reference genome of the specices you are working with