# exRNA.org data set search
it is very hard to use the data atlas. you can not do simple things like get sample counts, filter, group by, ...

- on data atlas page use the dataset tab

## Trending isolation kits
- 7,585 samples
- 7 named kits + "other"
-  80%  of sample
   * miRNeasy (1228) 
   * miRCury biofluids (Exiqon) (3265), 
   * MirVana PARIS (1539)
- we use Qiagen (603)

## Diverse Long RNAs Are Differentially Sorted into Extracellular Vesicles Secreted by Colorectal Cancer Cells

1. Diverse Long RNAs Are Differentially Sorted into Extracellular Vesicles Secreted by Colorectal Cancer Cells
   a. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6248336/ 11/21/2018
       i. previously showed that mutant KRAS colorectal cancer (CRC) cells release extracellular vesicles (EVs) containing distinct proteomes, microRNAs (miRNAs), and circular RNAs.
       ii. identify diverse classes of CRC, (colorrectal cancer),  extracellular long RNAs secreted in EVs and demonstrate differential export of specific RNAs
       iii. detected strong enrichment of Rab13 in mutant KRAS EVs and demonstrate functional delivery of Rab13 mRNA to recipient cells
       iv. implemented a CRISPR/ Cas9-based RNA-tracking system to monitor delivery to recipient cells
       v. show that gRNAs, (guide RNA) containing export signals from secreted RNAs can be transferred from donor to recipient cells. Our data support the existence of cellular mechanisms to selectively export diverse classes of RNA
       vi. TODO AEWIP check for intersting genes: Hinger et al. identify distinct subsets of cellular coding and long noncoding RNAs that are enriched in EVs that can be functionally transferred between cells, supporting a regulated form of cell-cell communication
       
   b. method
       i. utilized three isogenic colorectal cancer (CRC) cell lines that differ only in the mutational status of the KRAS gene (includes wild type)
       ii. previously showed that EVs from mutant KRAS CRC cells can be transferred to WT cells to induce cell growth, migration, and invasiveness
       iii. showed lncRNA can be functionaly transferred between cells
       
   c. results
       i.EVs derived from mutant KRAS cell lines have distinct proteomes, miRNA profiles, and circular RNA profiles compared to their parental cellular patterns and EVs from WT KRAS cells
       ii. Long RNA-seq was performed on rRNA-depleted total cellular RNA
           - majority of RNAseq reads were derived from rRNA
       iii. ??? performed correlation between replicants ???
       iv. analysis showed that specific RNAs are selectively exported into EVs
   d. data
       i. 18 cancer samples
       ii. [meta data](https://exrna-atlas.org/exat/gridview)
       iii RNA Isolation Kit: miRNeasy (Qiagen)
       iV human cell line
       V. biofluid: culture media, CONDITIONED
       vi. rna source: extracellular exosome
       vii. small RNA-seq
       vii. data links
        - [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA278673](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA278673)
        - [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67004](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67004]
        -
       
       
       
## Plasma extracellular RNA profiles in healthy and cancer patients
ref: [https://www.nature.com/articles/srep19413](https://www.nature.com/articles/srep19413)
- 2015


data
1. 192 samples 50 health, 142 cancer 
   - colorectal (50 male, 50 female) , prostate 36, pancreatic (3 male, 3 female)
2. biofluid: plasma
3. rna source: total cell-free biofluid RNA
4. RNA Isolation kit: miRNeasy (Qiagen)
5. profiling assay: small RNA-seq

- multivariate analysis of covariance, we identified significant associations
of these exRNAs with age, sex and different types of cancers
- developed multivariate statistical models to predict cancer status with an area under the
curve from 0.68 to 0.92 depending cancer type and staging

- def:  Exosomes are cell-derived membrane vesicles (30–100 nm)

-studies have shown that exRNAs are functionally active. For example, the release of extracellular
miRNAs is associated with anti-cancer signaling12 or cancer metastasis signaling13

- miR-143 within extracellular
vesicles inhibited proliferation of cancerous cells14 whereas introduction of miR-16 enriched extracellular vesicles
into prostate cancer cells significantly suppressed expression of miR-16 target genes15.


- issues (2015) with current studies
    * lack of systematic evaluation of the technical effect of RNA-seq library preparation on RNA abundance and
detectability.

    * lack of large scale RNA-seq data in a variety of populations including both healthy
individuals and those with diseases to generate exRNA expression profiles as a baseline reference.

- identified exRNA asscociated with age and set
  - is this something we need to control for?
- identified colorectal markers
- list most common mRNA found along with gene function
