

ref:
- [list of videos](https://statquest.org/video-index/)


# [intro to RNA-seq](https://www.youtube.com/watch?v=tlf6wYJrwKY&feature=youtu.be)

## step 1: prep rna-seq lib

1. isolate RNA
2. break into small fragments
   * illuminia can only deal wihtwith fragment lengths around 200 - 300 base pairs long
3. convert rna framents into double stranded DNA
   * dna is more stable than rna and can be easily amplified
4. add sequencing adaptors
    A. allow seq machine to recongomize framents
    b. allows you to seq different samples at the same time. different samples can use different adaptors. reduces cost
    c. does not work 100% of the time. ie. some sequence may not get adaptors
5. PCR amplification
   * only seq with adaptors are amplified; "they are enriched"
6. step 6 quality control
   a. verify library concentration
   b. verify library fragment lengths (not too long or short)
   
## step 2: sequence 
    * illiumina: flow cell, 'grid', each fragement is vertically oriented/attached on 'plate'flow cell plate
    * flurorescent probes nucleotide specific binding
    * machine takes picture
    * probe is washed off
    * probe bound to next necleotide in fragment
    

sequencing errors:

kinds of errors
- probe does not shine as bright as should
- lots of probs in same region with same color can lead to low quality score "low diversity" error
  * "low diversity error" are common when first few neculeoties are sequenced

machine generates quality score


Raw Data

- 4 lines for each rad
  * line 1: unique ID alway start with '@'
  * line 2: the bases called for fragment
  * line 3: always '+'
  * line 4: quality scores for each read
  
Next steps:

1. filter out garbage reads
   * low quality reads
   * reads that are artifacts of the chemistry
   * good read == adaptor seq + dna fragement + adapter seq
   * bad read == just adapater seqs
   
2. align the high quality reads to a reference genome
   * break genome into fragments
   * create index for genome fragments
   * for each read split into fragments and match to genome fragments
   * small reads make it easier to align when read is not not exactly same as ref genome
3. count the number of reads per gene
   * data file (bulk rna seq)
     + first col: gene name
     + remaining columns contain counts for each sample we sequenced
       - ie 3 cntrl and 3 treatment replicants
4. normalize read counts
   * need to account for difference between sample. i.e. one sample has more cells than another
   
   
## step 3Analyize data

1. plot the data
   * pca -> 2d, label cntr and treatment
     + are there outliers
     + expect to see biggest diff on x axis between cntrl & treatment
   * if sca, might want to exclused cells that are where seperation is not clear
2. identify differntially expressed genes between normal and mutant
   * [R], edgeR or DESeq2
   * some sort of volcan plot, vs logCPM vs logFoldChange
   * label point for interesting genes
   

