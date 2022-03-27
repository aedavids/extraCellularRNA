# Volcano Plots
quick hacks to create plots for lab research presenation on mond 3/28

1. copy data from mustard
   * used 1vsAllRunner.batch.sh to create data. See extracellularRNA/terra/deseq/R/README.md for instructions
   ```
   scp mustard:/scratch/aedavids/GTExData/data/deseq/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.te.gene.names.txt .
   
   scp mustard:/scratch/aedavids/GTExData/1vsAllRunner.sh/Validate_Lung_vs_all_results.csv .
   ```
   
2. start extraCellularRNA conda environment

3. add local python packages
   ```
   cd extraCellularRNA/terra/deseq/python
   export PYTHONPATH="${PYTHONPATH}:`pwd`/src"
   ```
   
4. run plots
   ```
   python plots/volcanoPlots.py -t "Validate Lung vs. not Lung n=119" \
       -i ../doc/plots/data/Validate_Lung_vs_all_results.csv \
       -o ../doc/plots/volcanoPlots/Validate_Lung_vs_all_results.png \
       -n 8 \
       -g ../doc/plots/data/gencode.v35.ucsc.rmsk.te.gene.names.txt
   ```

