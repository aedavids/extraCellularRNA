# Kras IPSC data overview

Follow up from separate meeting with Prof and Roman on 6/9/2020, and addition notes

Goal: find evidence that extracellular RNA influences IPSC differetiation

## data set :

## bulk intra cellular rna sequence sets

- fastq files (HiSeq read data)
    * /public/groups/kimlab/kras.ipsc/bulk.data/day.{5,7}/{ctrl,kras}.{1,2,3}/output_reverse_paired.fq.gz
    * /public/groups/kimlab/kras.ipsc/bulk.data/day.{5,7}/{ctrl,kras}.{1,2,3}/output_forward_paired.fq.gz
    * example first entry:
        + first line is read identifier
        + second line is sequence
        + third line is always '+'
        + forth line is read quality for each neculeotide
```
   (base) [aedavids@courtyard ctrl.1]$ zcat output_forward_paired.fq.gz | head -n 4
   
    @K00188:93:HVYFJBBXX:3:1101:1479:1490 1:N:0:GCCAAT
    GGCGTACAGCAGTGTCAACGTAGTAGTTAACAGGGTCTCCGCTGTGGATCATCAGGCCATCCACAAACTTCATGGATTTAGCCCTCTGTCCTCGGAGTTTC
    +
    AAFFFJAFJJJJJJJJJJ<-FFJJJJJJJJJJF7FJ-FFF7FJA-FFJFFJFFJJFJFJFFA7F7<A<JJJ<JJ<<FJFFFFJJJFFF7<<JJAAF<7A<J
```

- the alignement files
    * /public/groups/kimlab/kras.ipsc/bulk.data/day.{5,7}/bam.files/{ctrl,kras}.{1,2,3}star.out/pass.2/*.bam
    
- output files from star, salmon, can be found under
    * /public/groups/kimlab/kras.ipsc/bulk.data/day.{5,7}/{ctrl,kras}.{1,2,3}/
    
## bulk extra cellular data 'library created from media in petri dish'
we only have day 5 data

- fastq files (HiSeq read data)
    * /public/groups/kimlab/kras.ipsc/bulk.data/day.{5,7}/{ctrl,kras}.{1,2,3}/output_{forward,reverse}_{paried,unpaired}.fq.gz

- the output files form star, salmon, ... can be found under 
  * /public/groups/kimlab/kras.ipsc/exo.data/gen1c.day.5.exo.input.data/{ctrl,kras}.{1,2,3}
  
- the bam files
    * /public/groups/kimlab/kras.ipsc/exo.data/gen1c.day.5.exo.input.data/{ctrl,kras}.{1,2,3}/star.out/*.bam
    * bam example
      + first line
      
      ```
      samtools view ctrl.1.out.bam  | head -n 5 | tee ~/t
        $ cat ~/t | head -n 1 | tail -1
NB551675:38:HTMGLAFXY:1:11101:10183:1121    163 chr21   8254748 0   113M    =   8254748 113 GCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATT   AAAAAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEE<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE/EEEEEEEEEE6EEEEEEE<EAEEEEAEEAE   NH:i:8  HI:i:1  AS:i:224    nM:i:0
        ```
        
     + second line
        ```
        $ cat ~/t | head -n 2 | tail -1
        NB551675:38:HTMGLAFXY:1:11101:10183:1121    83  chr21   8254748 0   113M    =   8254748 -113    GCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATT   <///A//EAE///<A/EEEEEAE<<EEEEAEEEEEEEEEEEEAEEEEEEA/EE/EEEEEEEEEEAEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA   NH:i:8  HI:i:1  AS:i:224    nM:i:0
        
        
        ```


## intra cellular Single Cell
- /public/groups/kimlab/kras.ipsc/single.cell.data/day.5.ipsc.{ctrl,kras}
- /public/groups/kimlab/kras.ipsc/single.cell.data/day.7.eb.{ctrl,kras}
  * 'eb': embronic body, are three-dimensional aggregates of pluripotent stem cells

## extra cellular Single cell
- does not exist. Maybe not be possible to collect media from wells

Data sets we want to focus on as prep for meeting with Kevin are

# Daniel's directions:
a) The 3 IPSC scrambled exo rna libs
  - AED: only day 5, these are the control replicants

  - they are the extra cellular rna libs
  - sequenced rna from the dish not the cells
  - they are the control
  - the term 'scrambled' means
    * ref: p 105 BME notebook 9/30/19 - 'Srees examination prep lab meeting presentation'
    * IPSC cell culture is ? set up with cripsri system?
    * 20 nucliotide guide used to knock down kras
    * to ensure results are not due to added RNA as opposed to actual knock down
    * a random 20 nucliotied guide is used on the ctrl replicants. 
    * because it is random should not effect mRNA transcription
    
b) the 3 ipsc kras exo rna libs

    - they are extra cellular rna libs from the treatment. i.e. down regulated kras
    - sequenced rna from dish not the cells

Run desq to identify non biologic effect we need to control. E.G. any thing that is over or under expressed is probably due to experimental condition and not biology


## Make sure we know what data files to use
Roman mentioned that these where some of the first files he ran through his pipe line. 

"The only comparisons I am 110% on those we are submitting in the current manuscript. I have code that will run the in vitro comparisons no problem but I am currently writing a new script to handle the patient data and include metadata into the model. You can find my in vitro DESeq script @ /public/groups/kimlab/exoRNA-biomarkers-panc/scripts/tximport.and.normalize.R. R scripts are run with the ```Rscript``` command and behave like any other *nix script after that, it has argparse so you should be able to see what’s going on. "


```
(base) [aedavids@courtyard ~]$ echo $d
/public/groups/kimlab/kras.ipsc
```

The raw data files

```
(base) [aedavids@courtyard ~]$ ll $d/exo.data/bam.files
total 4.5G
-rw-rwxr-- 1 rreggiar kimlab 626M Feb 21 14:48 ctrl.1.out.bam*
-rw-rwxr-- 1 rreggiar kimlab 638M Feb 21 14:48 ctrl.2.out.bam*
-rw-rwxr-- 1 rreggiar kimlab 581M Feb 21 14:48 ctrl.3.out.bam*
-rw-rwxr-- 1 rreggiar kimlab 1.2G Feb 21 14:48 kras.1.out.bam*
-rw-rwxr-- 1 rreggiar kimlab 771M Feb 21 14:48 kras.2.out.bam*
-rw-rwxr-- 1 rreggiar kimlab 761M Feb 21 14:48 kras.3.out.bam*
-rw-rwxr-- 1 rreggiar kimlab 209K Feb 21 14:49 kras.alu.bed.gz*
```

### The "control libs

```
(base) [aedavids@courtyard ~]$ ls $d/exo.data/*ctrl*

# 1. the "type a" files
/public/groups/kimlab/kras.ipsc/exo.data/ctrl.de.seq.csv*


# 2. same data set as 1. used the desq normalized counts algo
/public/groups/kimlab/kras.ipsc/exo.data/ctrl.normalized.counts.csv*

# 3. same data as 1. just aligned to ???? reference includes transposible elements
/public/groups/kimlab/kras.ipsc/exo.data/ctrl.te.exo.comp.de.seq.csv*

# 4. same data as 3. used the desq normalized counts algo
/public/groups/kimlab/kras.ipsc/exo.data/ctrl.te.exo.comp.normalized.counts.csv*
```


### The treatment libs

```
(base) [aedavids@courtyard ~]$ ls $d/exo.data/*kras*

# 5. the "type b" files
/public/groups/kimlab/kras.ipsc/exo.data/kras.de.seq.csv*

# 6. same data as 5. used desq normalized count algo
/public/groups/kimlab/kras.ipsc/exo.data/kras.normalized.counts.csv*

# 7. same data as 5. aligned to ??? reference taht includes transposable elements
/public/groups/kimlab/kras.ipsc/exo.data/kras.te.exo.comp.de.seq.csv*

# 8. same as 7. used desq normalized counts algo
/public/groups/kimlab/kras.ipsc/exo.data/kras.te.exo.comp.normalized.counts.csv*
```


## Mastering DESQ2
Find a small set of files we can really dig into, improve the volcano plots
- label kras
- add color gradient to encode abundance info

Roman suggested using the 

```
public/groups/kimlab/kras.ipsc/day.5.de.seq.csv(base) [aedavids@courtyard ~]$ ll $d/day.{5,7}.de.seq.csv
# rna from the ipsc cells.
-rwxrwxr-x 1 rreggiar kimlab 2.2M Oct 30  2019 /public/groups/kimlab/kras.ipsc/day.5.de.seq.csv*
-rwxrwxr-x 1 rreggiar giuser 2.3M Oct 29  2019 /public/groups/kimlab/kras.ipsc/day.7.de.seq.csv*
```

avoid they do not look correct

```
(base) [aedavids@courtyard ~]$ ll $d/day.{5,7}.de.seq.norm.counts.csv
-rwxrwxr-x 1 rreggiar kimlab 3.6M Oct 30  2019 /public/groups/kimlab/kras.ipsc/day.5.de.seq.norm.counts.csv*
-rwxrwxr-x 1 rreggiar kimlab 3.7M Oct 29  2019 /public/groups/kimlab/kras.ipsc/day.7.de.seq.norm.counts.csv*
```
