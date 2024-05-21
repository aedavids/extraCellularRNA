# Lump Brain samples
andrew E. Davidson  
aedavids@ucsc.edu
5/21/24

In general the 13 brain types have poor sensitivity. Try lumping them. The advantage of lumping is
1. improve brain classification 
2. reduce the number of signature genes
3. improve over all classificaiton rates
4. should be easier than dropping the brain samples. We only need to run 1 vs. all for the new lump


## Step 1. run "brain" 1 vs. all

a. create new colData. We can use sed
    ```
    cd /private/groups/kimlab/GTEx_TCGA/groupbyGeneTrainingSets
    cat GTEx_TCGA_TrainColData.csv | sed  's/\([^,]*,[^,]*,\)Brain_[^,]*/\1Brain/' > GTEx_TCGA_TrainLumpBrain.colData.csv
    ```

    check
    ```
    $ grep Brain GTEx_TCGA_TrainColData.csv | cut -d , -f 3 | sort | uniq -c 
     91 Brain_Amygdala
    106 Brain_Anterior_cingulate_cortex_BA24
    148 Brain_Caudate_basal_ganglia
    129 Brain_Cerebellar_Hemisphere
    144 Brain_Cerebellum
    153 Brain_Cortex
    126 Brain_Frontal_Cortex_BA9
    118 Brain_Hippocampus
    121 Brain_Hypothalamus
    147 Brain_Nucleus_accumbens_basal_ganglia
    123 Brain_Putamen_basal_ganglia
     95 Brain_Spinal_cord_cervical_c-1
     84 Brain_Substantia_nigra
     
    $ grep Brain GTEx_TCGA_TrainLumpBrain.colData.csv | cut -d , -f 3 | sort | uniq -c 
   1585 Brain
    ```
b. create a new 1vsAllLumpBrain dir
    - symbolically link  /private/groups/kimlab/GTEx_TCGA/1vsAll/* 
    - delete results startingn with "Brain"
    
        ```
        $ cd /private/groups/kimlab/GTEx_TCGA
        $ mkdir 1vsAllLumpBrain
        $ cd 1vsAllLumpBrain
    
        $ ln -s ../1vsAll/*.results .
        
        $ rm Brain*
        ```
    
c. create a driver scripts
    - ref : extraCellularRNA/intraExtraRNA_POC/adenocarcinoma.vs.control/run.adenocarcinoma.vs.control.sh
      ```
      $ d=/private/home/aedavids/extraCellularRNA/intraExtraRNA_POC/adenocarcinoma.vs.control
      $ cp $d/adenocarcinoma.vs.control.1vsAllTask.input.json brain.vs.all.input.json
      $ cp /private/home/aedavids/extraCellularRNA/intraExtraRNA_POC/adenocarcinoma.vs.control/adenocarcinoma.cromwellOptions.json .
      $  mv adenocarcinoma.cromwellOptions.json brain.cromwellOption.json
      $ cp $d/adenocarcinoma.vs.control/run.adenocarcinoma.vs.control.sh run.brain.vs.all.sh
      ```
      
   - edit
     *  brain.cromwellOption.json: 
        _ change final_workflow_outputs_dir
    *  brain.vs.all.input.json
       - deseq_one_vs_all.one_vs_all.referenceLevel : Brain
       - deseq_one_vs_all.one_vs_all.colData
         + /private/groups/kimlab/GTEx_TCGA/groupbyGeneTrainingSets/GTEx_TCGA_TrainLumpBrain.colData.csv
       - deseq_one_vs_all.one_vs_all.design
         + "~  sex + tissue_id"
       - deseq_one_vs_all.one_vs_all.countMatrix
         + "/private/groups/kimlab/GTEx/GTExTrainGroupByGenesCountMatrix.csv"
   * run.brain.vs.all.sh
     - --inputs
       +  brain.vs.all.input.json
     - --options
       + brain.cromwellOption.json
     
 d. execute run.brain.vs.all.sh
   ```
   cd /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/lumpBrain
   $ setsid sh -c 'set -x;run.brain.vs.all.sh' > run.brain.vs.all.sh.out 2>&1 &
   ```
    
## Step 2. run  deconvolution hyper parameter tunning pipeline model best500FindAllDegree1_wl500

- this will crash and burn. we only need the intersection dictionary created by the upsetPlot code

## Step 3. run deconvolution hyper parameter tunning pipeline model after best10CuratedDegree1.s

- create a new driver script. Model after extraCellularRNA/deconvolutionAnalysis/bin/1vsAll-~gender_category/best10CuratedDegree1.sh. This run is one of our best results. It also automatically select the degree 1 genes to use. It does not rely on curatedGeneSignature

- change colData
  * use  our lumped colData

- change deseqResultsDir
  * use 1vsAllLumpBrain dir
  * upstream pipeline reads all the deseq results file in the directory and parses the file names to determin which tissue and cancer type are present
  
- change upstreamRun use our lumped brain best500FindAllDegree1_wl500

- change curatedGeneSignature use our lumped brain best500FindAllDegree1_wl500
