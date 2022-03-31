# Data overview 3/28/22
At this point in time we have the numReads matix and the groupby matrix for GTEx.
We want to do a priliminary exploreation about our gene signature method by looking
at a few tissue types

```
ls /private/groups/kimlab/GTEx/
GTExTestColData.csv                   GTExTrainNumReadsMatrix.tsv
GTExTestGroupByGenesCountMatrix.csv   GTExValidateColData.csv
GTExTestNumReadsMatrix.tsv            GTExValidateGroupByGenesCountMatrix.csv
GTExTrainColData.csv                  GTExValidateNumReadsMatrix.tsv
GTExTrainGroupByGenesCountMatrix.csv
```

- we have a copy of these files on mustart /scratch/aedavids/GTExData. 
- We ran startRStudioServer.sh to start up our 1vsAll docker container.
- we edited and ran 1vsAllRunner.batch.sh to create some preliminary data.
  - used the validate data sets. it is much smaller than training. it was faster
  ```
  docker exec --detach --user rstudio  adoring_darwin \
      /home/rstudio/extraCellularRNA/terra/deseq/R/1vsAllRunner.batch.sh
  ```

Select data sets where copied back to the data dir. It is easier to create plots on mac.

# Volcano Plots
quick hacks to create plots for lab research presenation on mond 3/28

1. copy data from mustard
   * used 1vsAllRunner.batch.sh to create data. See extracellularRNA/terra/deseq/R/README.md for instructions
   ```
   scp mustard:/scratch/aedavids/GTExData/data/deseq/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.te.gene.names.txt .
   
   scp mustard:/scratch/aedavids/GTExData/1vsAllRunner.sh.out/Validate_Lung_vs_all_results.csv .
   ```
   
2. start extraCellularRNA conda environment

3. add local python packages
   ```
   cd extraCellularRNA/terra/deseq/python
   export PYTHONPATH="${PYTHONPATH}:`pwd`/src"
   ```
   
4. run plots easy way to run many plots by just chaning refLevel and n variables
   ```
   refLevel=Thyroid; 
   n=130; 
   dataRoot=../doc/plots/data; 
   inputFile="${dataRoot}/Validate_${refLevel}_vs_all_results.csv"; 
   outputFile="${dataRoot}/../volcanoPlots/Validate_${refLevel}_vs_all_results.png"; 
   title="Validate ${refLevel} vs. not ${refLevel} n=$n";
    
    python plots/volcanoPlots.py \
        -n $n \
        -i $inputFile \
        -o $outputFile \
        -g ../doc/plots/data/gencode.v35.ucsc.rmsk.te.gene.names.txt \
        -t "$title"
   ```

# Upset Plots
shows intersection between tissue signature profiles

1. copy the 1vsAll results from mustard
   ```
   scp -r mustard:/scratch/aedavids/GTExData/1vsAllRunner.sh.out/ .
   ```

2. on mac edit and run 
   - extraCellularRNA/terra/deseq/doc/plots/jupyterNotebooks/GTExValidateExploration.ipynb.
   - it will select the ?best? 25 signature genes for each tissue type


3. create dataSet.csv file 
   ```
   cd extraCellularRNA/terra/deseq/doc/plots/data/1vsAllRunner.sh.out/design:~+sex+tissueId
   
   createDataSetCSV.sh
   ```
   

4. start extraCellularRNA conda environment and add upset plot package
   ```
   cd ~andrewdavidson/googleUCSC/kimLab/unmappedReadsAnalysis/python
   export PYTHONPATH="${PYTHONPATH}:`pwd`"
   ```
   
5. create plots
   ```
   cd extraCellularRNA/terra/deseq/python
   python plots/geneSignatureUpsetPlot.py \
       -t "my title" \
       -n 8  \
       -i ../doc/plots/data/signatureGenes*.csv \
       -o ../doc/plots/upsetPlots/signatureGenes.png
       
   python plots/geneSignatureUpsetPlot.py \
       -t "Validation Set: n=25 Signature Genes, padj < 0.001, lf2c > 2, sorted by baseMean" \
       -d ../doc/plots/data/1vsAllRunner.sh.out/design_~+sex+tissueId/dataSets.csv \
       -o ../doc/plots/foo.upsetPlots/signatureGenes.png 
   ```
   
# bar charts

```
python plots/barChart.py -t "validate set"  \
    -i ../doc/plots/data/validateTissueIds.csv \
    -o ../doc/plots/barCharts/validateTissueTypes.png
```

