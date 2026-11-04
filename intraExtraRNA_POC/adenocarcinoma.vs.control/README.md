# run 1vsAll on nanopore Adenocarcinoma samples
aedavids@ucsc.edu  
5/4/24  


**overview**  
our random forest model trained on elife ["Healthy donor", "Esophagus Cancer"] samples and ESCA biomarkers does not preform well on our nanopore Adenocarcinoma. **The problem may be the control samples are from really old people. They where intended to be used with alzheimer's not Adenocarcinoma**

extraCellularRNA/intraExtraRNA_POC/jupyterNotebooks/elife/elifeBinaryRandomForestResults.ipynb

extraCellularRNA/londonCalling2024/jupyterNotebooks/nanoporeAdenocarcinomaBinaryClassification.ipynb

HUGO_Genes
['AC012615.3', '(TA)n', 'UBE2SP2', 'HERVFH19-int', 'PRELID1P1', 'LTR106', 'AC010336.3', 'GOLGA8S', 'MER5C', 'CCDC160']

elifeLungGenes
['LTR106_Mam', 'ENSG00000217325.2', 'ENSG00000224126.2', 'ENSG00000203952.9', 'ENSG00000267125.2', 'MER5C', 'HERVFH19-int', 'ENSG00000268120.1', '(TA)n', 'ENSG00000261739.2']

## try adding additional biomarkers

### 1)  1vsAll on nanopore Adenocarcinoma samples
ref: extraCellularRNA/terra/natureBioMedEng  

a) create a working director and copy a template WDL workflow input json file
```
cd extraCellularRNA/intraExtraRNA_POC

mkdir adenocarcinoma.vs.control

rootDir=/private/home/aedavids/extraCellularRNA/terra/natureBioMedEng
cp ${rootDir}/GTExWhole_Blood.vs.GTEx/GTExWhole_Blood.vs.GTEx.1vsAllTask.input.json adenocarcinoma.vs.control.1vsAllTask.input.json
```

b) edit the cromwell input config file  adenocarcinoma.vs.control.1vsAllTask.input.json.

The count file only contains control and Adenocarcinoma samples. 

```
{
    "deseq_one_vs_all.one_vs_all.memoryGb": "64",
    "deseq_one_vs_all.one_vs_all.runTimeCpu": "2",
    "deseq_one_vs_all.one_vs_all.colData": "/private/groups/kimlab/aedavids/londonCalling2024/data/colData.high.adenocarcinoma.control.csv", 
    "deseq_one_vs_all.one_vs_all.isCSV": "true",
    "deseq_one_vs_all.one_vs_all.isDebug": "true",
    "deseq_one_vs_all.one_vs_all.referenceLevel": "adenocarcinoma",
    "deseq_one_vs_all.one_vs_all.diskSpaceGb": "80",
    "deseq_one_vs_all.one_vs_all.design": "~  sex + histology",
    "deseq_one_vs_all.one_vs_all.dockerImg": "aedavids/edu_ucsc_kim_lab-1vsall_1.0",
    "deseq_one_vs_all.one_vs_all.runTimePreemptible": "1",
    "deseq_one_vs_all.one_vs_all.countMatrix": "/private/groups/kimlab/aedavids/londonCalling2024/data/normalizedCounts.high.adenocarcinoma.control.csv"
}

```

c) edit cromwellOptions.json config file adenocarcinoma.cromwellOptions.json . It controls where output will be written
log files will be written to current working directory. ie the dir you ran the run script from 

```
{
    "final_workflow_outputs_dir": "/private/groups/kimlab/aedavids/londonCalling2024/run.adenocarcinoma.vs.control.sh.out",    
    "use_relative_output_paths": true
}
```

d) run using setsid

```
setsid sh -c 'set -x;run.adenocarcinoma.vs.control.sh' > run.adenocarcinoma.vs.control.sh.out 2>&1 &

```

see ~/extraCellularRNA/terra/natureBioMedEng/README.md for how to monitor, debug, ...

### 2) enrich our deconvolution biomarkers. 
We currently use degree 1 genes found in window length of 500. consider adding ESCA degree 2 genes. These genes are less likely to add ambiguity
our best hyperparmeter tunning results are /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best10CuratedDegree1_ce467ff created by best10CuratedDegree1.sh it searched the upset plot intersection dictionary /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500FindAllDegree1_wl500/training/best500FindAllDegree1_wl500.sh.out/upsetPlot.out/best500_findAllDegree1_wl500.intersection.dict



### 3) explore best100 ESCA volcano plots
see /private/home/aedavids/extraCellularRNA/londonCalling2024/jupyterNotebooks/ESCA1vsAllVolcanoPlots.out


