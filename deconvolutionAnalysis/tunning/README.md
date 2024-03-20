# Gene Signatue Profile fine tunning

Andrew E. Davidson  
<aedavids@ucsc.edu>  
1/1/24

best100Enriched_6_Degree1GTEx_TCGA LUAD, LUSC sensitivity is 0.673 and 0.704.  When we look at the false positives and false negatives for these two classes it looks like we can 31 miss classifications. Reducing these error should improve sensitivity.

run LUAD vs. LUSC to find good discriminators.  

ref: extraCellularRNA/terra/natureBioMedEng  

## install

```
rootDir=/private/groups/kimlab/aedavids/deconvolution/LUAD.vs.LUSC
mkdir -p $rootDir

cd extraCellularRNA/deconvolutionAnalysis/tunning

files="cromwellOptions.json LUAD.vs.LUSC.1vsAllTask.input.json LUSC_LUAD_ColData.csv LUSC_LUAD_GroupByGenseCountMatrixData.csv runLUAD.vs.LUSC.1vsAll.sh"

cp $files $rootDir
```

## run

run in new session and collect std out and std error. You can logout. script will continue to run

```
cd $rootDir
logFile="${rootDir}/runLUAD.vs.LUSC.1vsAll.sh.log"
setsid sh -c 'set -x; runLUAD.vs.LUSC.1vsAll.sh' > ${logFile} 2>&1 & 
```
