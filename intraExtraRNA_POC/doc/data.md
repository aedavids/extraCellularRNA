# Intracellular, Extracellular RNA proof of concept Data overview
Andrew E. Davidson
aedavids@ucsc.edu


**Data Sets**
- [Roman's normalized data sets](#roman)
- [kimlab elif counts](#elife)
- [panc.plasma.2020](#panc.plasma.2020)
- [pancreas.plasma.ev.long.RNA](#pancreas.plasma.ev.long.RNA)
- [<span style="color:red"> what is panc_resub?</span>](#panc_resub)

<hr />
<a name="roman">**Roman's normalize data sets**</a>

**elife**
```
cd /private/groups/kimlab/aedavids_data
ls elife *
elife_all_norm_counts_2023-05-18.csv
elife_scaled_metaData_2023-05-18.csv

aedavids@mustard $ wc -l elife*
    76556 elife_all_norm_counts_2023-05-18.csv
      225 elife_scaled_metaData_2023-05-18.csv
    76781 total
```

meta data columns
```
aedavids@mustard $ head -n 1 elife_scaled_metaData_2023-05-18.csv 
"","salmon_frag_length_mean","salmon_frag_length_sd","salmon_num_processed",
"salmon_num_mapped","salmon_num_decoy_fragments","salmon_num_dovetail_fragments",
"salmon_num_fragments_filtered_vm",
"salmon_num_alignments_below_threshold_for_mapped_fragments_vm","salmon_percent_mapped",
"age","diagnosis","gender","input_vol","dataset","condition"
```

example meta data
```
$ head -n 3 elife_scaled_metaData_2023-05-18.csv | cut -d , -f 1-5
"","salmon_frag_length_mean","salmon_frag_length_sd","salmon_num_processed","salmon_num_mapped"
"SRR14506659",-0.74719344483705,-0.575976877708043,-0.842040746773851,-1.42249054813429
"SRR14506660",0.260123326819362,0.410348155985916,-0.701701964699133,-0.865153222330587
```

```
$ head -n 3 elife_scaled_metaData_2023-05-18.csv | cut -d , -f 6-9
"salmon_num_decoy_fragments","salmon_num_dovetail_fragments","salmon_num_fragments_filtered_vm","salmon_num_alignments_below_threshold_for_mapped_fragments_vm"
-1.24676754658216,-0.945323019739035,-1.0770993514578,-0.429068029826966
0.0099671790600485,-0.173367717315595,-0.141121784435817,-0.160291452098426
```

```
$ head -n 3 elife_scaled_metaData_2023-05-18.csv | cut -d , -f 10-16
"salmon_percent_mapped","age","diagnosis","gender","input_vol","dataset","condition"
-2.00444714324715,-0.258071341194335,"Esophagus Cancer","male",NA,"elife","Esophagus Cancer"
-0.779962818605958,0.519178815814486,"Esophagus Cancer","male",NA,"elife","Esophagus Cancer"
```

224 samples
```
aedavids@mustard $ head -n 1 elife_all_norm_counts_2023-05-18.csv | comma2newLine | wc -l
224

aedavids@mustard $ head -n 1 elife_all_norm_counts_2023-05-18.csv | comma2newLine | head -n 3
SRR14506659
SRR14506660
SRR14506661
```

<hr />

**panc**

```
ls panc_normal_*
panc_normal_norm_counts_2023-05-18.csv  panc_normal_scaled_metaData_2023-05-18.csv
aedavids@mustard $ 

$ wc -l panc_normal_*
   76556 panc_normal_norm_counts_2023-05-18.csv
      33 panc_normal_scaled_metaData_2023-05-18.csv
```

meta data columns
```
$ head -n 1 panc_normal_scaled_metaData_2023-05-18.csv 
"","sample","salmon_frag_length_mean","salmon_frag_length_sd","salmon_num_processed","salmon_num_mapped","salmon_num_decoy_fragments","salmon_num_dovetail_fragments","salmon_num_fragments_filtered_vm","salmon_num_alignments_below_threshold_for_mapped_fragments_vm","salmon_percent_mapped","gender","age","diagnosis","stage","dataset","input_vol","condition"
```

example meta data
```
$ head -n 3 panc_normal_scaled_metaData_2023-05-18.csv | cut -d , -f 1-5
"","sample","salmon_frag_length_mean","salmon_frag_length_sd","salmon_num_processed"
"1","panc.1.2.3",2.21748201920628,2.36138106999858,-1.18754005718277
"2","panc.2.2.7",-0.447018603746171,0.584389506497176,-1.22336768997631
aedavids@mustard $ 

$ head -n 3 panc_normal_scaled_metaData_2023-05-18.csv | cut -d , -f 6-8
"salmon_num_mapped","salmon_num_decoy_fragments","salmon_num_dovetail_fragments"
-1.06505555478815,-0.993900533145792,-1.55168941465303
-1.11422968803571,-0.864212601544075,-1.60689306047192
aedavids@mustard $ 

$ head -n 3 panc_normal_scaled_metaData_2023-05-18.csv | cut -d , -f 9-12
"salmon_num_fragments_filtered_vm","salmon_num_alignments_below_threshold_for_mapped_fragments_vm","salmon_percent_mapped","gender"
-1.30336783091096,-0.391112172078636,1.02396636772035,"Male"
-1.23780982673013,0.109589116721385,0.453104655498845,"Male

$ head -n 3 panc_normal_scaled_metaData_2023-05-18.csv | cut -d , -f 13-18
"age","diagnosis","stage","dataset","input_vol","condition"
-0.167944278512434,"panc","IV","pilot",1.21949365854475,"panc"
2.36109897438068,"panc","IV","pilot",1.5621834242829,"panc"
```

num samples 32
```
$ head -n 1 panc_normal_norm_counts_2023-05-18.csv | comma2newLine | wc -l
32

 $ head -n 1 panc_normal_norm_counts_2023-05-18.csv | comma2newLine | head -n 3
panc.1.2.3
panc.2.2.7
panc.3.2.5

$ head -n 1 panc_normal_norm_counts_2023-05-18.csv | comma2newLine | tail  -n 3
n43
n44
n45
```

<hr />
<a name="elife">**kim lab elife data set**</a>

- see p 129 5/9/2023 "roman intra/extra RNA POC" in notebook jan 14 2022
- Roman has normlized complete seq data for all of the elif samples.
- control + colorectal cancer, stomach cancer, liver cancer, lung cancer, and esophageal cancer

<hr />
<a name="panc.plasma.2020">**/private/groups/kimlab/panc.plasma.2020/data**</a>

17 control samples 14 panc

Do we know where these sample came from?

Do we have any patient data data?

what is the total size of the data set?

The salmon cmd_info.json contains the command line arguments
```
cd /private/groups/kimlab/panc.plasma.2020/data

$ ls
all.plasma.data/  unprocessed.plasma.data/
```

```
ls all.plasma.data/
ctrl.10.0.4/  ctrl.1.3.0/   ctrl.5.1.3/   panc.10.2.5/  panc.6.1.1/  panc.9.2.5/
ctrl.1.0.5/   ctrl.13.0.7/  ctrl.6.1.5/   panc.1.2.3/   panc.6.2.5/
ctrl.11.1.1/  ctrl.1.5.0/   ctrl.7.1.2/   panc.2.2.7/   panc.7.1.1/
ctrl.1.1.5/   ctrl.2.3.4/   ctrl.8.1.2/   panc.3.2.5/   panc.7.2.5/
ctrl.1.2.0/   ctrl.3.2.5/   ctrl.9.1.2/   panc.4.3.9/   panc.8.2.5/
ctrl.12.1.1/  ctrl.4.1.3/   panc.10.1.2/  panc.5.4.5/   panc.9.1.1/

$ ll all.plasma.data/panc.10.2.5/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/quant.sf
-rw-r--r-- 1 aedavids prismuser 606M Jun 18  2021 all.plasma.data/panc.10.2.5/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/quant.sf
```

```
$ ls unprocessed.plasma.data/
pah.118.vol/  pah.220.vol/  pah.423.vol/  pah.49.vol/  pah.577.vol/
pah.154.vol/  pah.222.vol/  pah.441.vol/  pah.52.vol/  pah.72.vol/

$ ll unprocessed.plasma.data/pah.118.vol/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/quant.sf
-rw-r--r-- 1 aedavids prismuser 605M Jun 18  2021 unprocessed.plasma.data/pah.118.vol/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx/quant.sf
```

<hr/>
<a name="pancreas.plasma.ev.long.RNA">**/private/groups/kimlab/pancreas.plasma.ev.long.RNA**</a>

- 'Plasma extracellular vesicle long RNA profiling identifies a diagnostic signature for the detection of pancreatic ductal adenocarcinoma'
- [paper](https://gut.bmj.com/content/69/3/540)
- README.md . contains citation, location of meta data, ...

```
ls data/healthy/ | wc -l
117

ls data/PDAC/ | wc -l
284

data/PDAC/SRR10080544/salmon.out/quant.sf
data/PDAC/SRR10080544/SRR10080544_pass_1.fastq.gz
data/PDAC/SRR10080544/SRR10080544_pass_2.fastq.gz
```

**<span style="color:red">what? how? want normalization done. What is diff between two dirs?</span>**
```
ls -1 data/normalizedCounts/normalizedCounts.gencode.v35.tx.to.gene.csv/
pancreas.plasma.ev.long.RNA.normalized.deseq.BIO_TYPE.counts.csv
pancreas.plasma.ev.long.RNA.normalized.deseq.ENST.counts.csv
pancreas.plasma.ev.long.RNA.normalized.deseq.HGNC.counts.csv
```

```
ls -1 data/normalizedCounts/normalizedCounts.gencode.v35.ucsc.rmsk.tx.to.gene.csv/
pancreas.plasma.ev.long.RNA.normalized.deseq.biotype.counts.csv
pancreas.plasma.ev.long.RNA.normalized.deseq.gene.counts.csv
pancreas.plasma.ev.long.RNA.normalized.deseq.tx.counts.csv
```

<hr />
<a name="panc_resub">**/private/groups/kimlab/panc_resub**</a>

<span style="color:red">?????? </span>

