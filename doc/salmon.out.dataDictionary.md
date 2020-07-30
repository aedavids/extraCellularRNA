# Salmon output data dictionary

sample output

```
gi = "/public/groups/kimlab/kras.ipsc/bulk.data/day.7/ctrl.1/gencode.salmon.out/"
dataRoot = pl.Path(gi)
giData = 'quant.sf'
salmonFilePath = dataRoot.joinpath(giData)
salmonDF = pd.read_csv(salmonFilePath, delimiter='\t')
salmonDF.head()
```

```
Name	Length	EffectiveLength	TPM	NumReads
0	ENST00000456328.2|ENSG00000223972.5|OTTHUMG000...	1657	674.715	0.422396	11.070
1	ENST00000450305.2|ENSG00000223972.5|OTTHUMG000...	632	451.000	0.000000	0.000
2	ENST00000488147.1|ENSG00000227232.5|OTTHUMG000...	1351	988.761	28.248736	1084.876
3	ENST00000619216.1|ENSG00000278267.1|-|-|MIR685...	68	8.000	0.000000	0.000
4	ENST00000473358.1|ENSG00000243485.5|OTTHUMG000...	712	519.163	0.212906	4.293
```

## Understanding the 'name' column
two examples of "name" col values from a salmon output file

```
ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|lncRNA|

ENST00000450305.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000002844.2|DDX11L1-201|DDX11L1|632|transcribed_unprocessed_pseudogene|
```

- There are multiple name delimintated by '|'.
- col 0: ENST00000450305.5
    * Ensembl transcript stable id with version information. 
    * the version information is '.5'
- col 1: ENSG00000223972.5
    * Ensembl gene stable id with version information
    * the version information is '.5'
- col 2: OTTHUMG00000000961.2
    * ???? gene name ??? with version information
- col 3: OTTHUMT00000362751.1
    * ???? transcript name with version information ???
- col 4: DDX11L1-202
    * DDX11L1 is HGNC Symbol '202' is a specific transcript for the gene. Some time the term iso form is used instead of transcript
- col 5: DDX11L1
    * HGNC Symbol
- col 6: 1657
    * length of reference transcript. do not confuse this with transcript length col.
- col 7: reference bio type
    
bio type counts
```
# n = 1 return first split
# n = -1 return all splits
# expand = True : return DataFrame/MultiIndex expanding dimensionality
pd.set_option('max_colwidth', None) # no limit
newNamesDF = salmonDF["Name"].str.split("|", n = -1, expand = True) 
newNamesDF.iloc[:, 7].value_counts()
```
 
```
protein_coding                        83728
lncRNA                                74961
retained_intron                       28383
nonsense_mediated_decay               15788
processed_pseudogene                  10143
unprocessed_pseudogene                 2596
misc_RNA                               2175
snRNA                                  1836
miRNA                                  1828
TEC                                    1145
snoRNA                                  933
transcribed_unprocessed_pseudogene      926
transcribed_processed_pseudogene        493
rRNA_pseudogene                         488
IG_V_pseudogene                         185
IG_V_gene                               143
transcribed_unitary_pseudogene          138
TR_V_gene                               106
unitary_pseudogene                       98
non_stop_decay                           91
TR_J_gene                                78
polymorphic_pseudogene                   63
scaRNA                                   48
pseudogene                               38
TR_V_pseudogene                          33
IG_D_gene                                30
rRNA                                     23
IG_C_gene                                23
Mt_tRNA                                  22
IG_J_gene                                18
IG_C_pseudogene                           9
ribozyme                                  8
TR_C_gene                                 6
sRNA                                      5
TR_D_gene                                 4
TR_J_pseudogene                           4
IG_J_pseudogene                           3
Mt_rRNA                                   2
translated_processed_pseudogene           2
translated_unprocessed_pseudogene         2
vaultRNA                                  1
scRNA                                     1
IG_pseudogene                             1
```
