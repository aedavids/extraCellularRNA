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

