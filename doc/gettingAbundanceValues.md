# Get getting Abundance values
Goal: color volcano plots to identify signifigant over/under expressed genes with high read counts. 

initial data to focus on. Goal really understand dseq2

## How where the original DESQ2 files created?

These files are kras ipsc ?? bulk, not SCA ??? 
```
kras.ipsc.day.5.de.seq.csv 
kras.ipsc.day.7.de.seq.csv
```

There are three kras libraries
```
(base) [aedavids@courtyard day.7]$ ls $d/
bam.files/  ctrl.1/  ctrl.2/  ctrl.3/  kras.1/  kras.2/  kras.3/
```

<span style="color:red"> What does the desq data mean? Deseq is a ratio? </span> 
There are three replicaints/libries yet only one set of metrics for each gene. 

```
(extraCellularRNA) $ head kras.ipsc.day.7.de.seq.csv
"","baseMean","log2FoldChange","lfcSE","pvalue","padj"
"5S",3.03435747948001,0.000507678899665888,0.329042988338836,0.990510389964409,NA
"A1BG",136.214097908835,-0.0082954500999119,0.189893416836959,0.956030632440566,0.980079512557846
```

Roman mentioned we can not use the  normalized counts produced by desq2. There originally run was buggy

## Lab note from Roman
- master pipeline scripts 
  /public/groups/kimlab/exoRNA-biomarkers panc/scripts/patient.tximport.and.normalize.R
  
- count data 
  The abundance values for initial import into DESEQ? Those are the quant.sf files in the salmon output directories of interest
  ```
  find . -type d -name "*salmon*" -print | tee ~/findSalmon.out
  ```
  
  ```
  (base) [aedavids@courtyard gencode.salmon.out]$ ls -l $d/quant.sf
-rw-rwxr-- 1 rreggiar giuser 604M Feb 29 02:13 /public/groups/kimlab/kras.ipsc/bulk.data/day.7/kras.1/gencode.salmon.out/quant.sf*
  ```
  
  <span style="color:red"> looks like we need to convert to hugo names? Can use TPM or do we need to normalize some how? The values look strange. fractional TPM looks okay but how can numReads be fractional?</span>
  
  4:00 Roman, "the name is gene code it has all possible names, we should find the hugo name in there" 
  
  "salmon does kmer from reads to kmers of transcriptome. hence the numReads are not integer"

  ```
  (base) [aedavids@courtyard gencode.salmon.out]$ head quant.sf 
  Name	Length	EffectiveLength	TPM	NumReads
  ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|lncRNA|	1657	696.686	0.217458	2.765
  ```
  
# Notes from vignette
ref: beginers guide to desq2 http://biorxiv.org/lookup/doi/10.1101/002832

## preparing count matrix
The value in the i-th row and the j-th column of the matrix tells how many reads have been mapped to gene i in sample j.

2.1 counts must be raw counts
2.2 align reads
    - do we need to sort?
2.3 count reads
    - need to konw if strand specific or not
    - what package? program did we use?
    - where is the pipe line script ?
    - where is the output?
    - where are the gene model gtf file?
    ```
    /public/groups/kimlab
    find . -name "*.gtf" -print | tee ~/findGTF.out
    base) [aedavids@courtyard kimlab]$ grep -v Perm ~/findGTF.out
    ./aale.kras/data/other/wetlab/extras/gencode.v32.annotation.gtf
    ./aale.kras/data/other/wetlab/extras/sex.linked.anno.gtf
    ./kras.ipsc/nanopore.data/flair.collapse.isoforms.gtf
    ./genomes.annotations/gencode.te.cons.gtf
    ./genomes.annotations/GRCh38_rmsk_TE.gtf
    ./genomes.annotations/ltr7.gtf
    ./genomes.annotations/gencode.31/gencode.v31.annotation.gtf
    ./genomes.annotations/gencode.32/gencode.v32.annotation.gtf
    ```
    - is the "sumurized data stored on disk?
    
    
### Can we find the raw count data input to Desq2?

looks like the reads where alligned using salmon.alevin.sh, 

```
(base) [aedavids@courtyard kras.ipsc]$ pwd
/public/groups/kimlab/kras.ipsc
(base) [aedavids@courtyard kras.ipsc]$ !find
find -type d  -name alevin.out -print
./single.cell.data/DKSL010A/alevin.out
./single.cell.data/DKSL010B/alevin.out
./single.cell.data/DKSL010C/alevin.out
./single.cell.data/DKSL010D/alevin.out
(base) [aedavids@courtyard kras.ipsc]$
```

???
```
(base) [aedavids@courtyard alevin]$ pwd
/public/groups/kimlab/kras.ipsc/single.cell.data/DKSL010B/alevin.out/alevin
(base) [aedavids@courtyard alevin]$ zcat quants_mat.mtx.gz| head -n 10
%%MatrixMarket	matrix	coordinate	real	general
2407	59033	5007650
1	2	0.333333
1	33	1.170588
1	36	0.863179
1	38	19.745705
1	59	3.000000
1	68	2.000000
1	84	2.000000
1	88	7.000000
```

next step is to sort

then count 

```
(base) $ head bad.day.5.de.seq.norm.counts.csv 
"","ctrl.1","ctrl.2","ctrl.3","kras.1","kras.2","kras.3"
"5_8S_rRNA",0,0,0,0,0,0
"5S",0,0,0,0,0,0
"5S_rRNA",0,0,0,0,0,1.13219185921882
"7SK",0,0,0,0,1.18644898909649,0
"7SLRNA",0,0,0,0.936152362010547,0,0
"A1BG",153.6231373634,106.424870046778,110.083539326623,114.984081937783,132.445131827696,132.776017879379
"A1BG-AS1",86.0372735133266,91.839408506211,85.8369106338823,112.434887327076,96.2609115030601,92.8753998547189
"A1CF",0,0,0,0,0,0
"A2M",5.36803798907775,1.80417656899995,0.925520474233922,5.04222502242178,0,4.02811450670747
(base) $ 

```
