# Clean Up

We need to reduce our data storage costs

## Uber, can we delete
- data/matrices/NumReads/TCGA_UVM_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab.NumReads.tsv
  * 33 tcga tsv files that where created by https://app.terra.bio/#workspaces/test-aedavids-proj/uber/analysis/launch/createCountMatrix.ipynb
  


## Clean up fail TCGA salmon quant runs 
- a) generateTCGAMatrixCreationScripts.ipynb 4/27/22
  * run on mac
  * creates script used to transfer salmon quant files to jupyter notebook vm running in terra
  * finds results from 'salmonTarQuantWorkflow v 4'
    + TCGA workspace data varies from cohort/study group to studyGroup
    + bam files may be in tar, or gzip, contain paired or single end reads, and include replicants
    + It was tricky and expensive to sort this out. for each data model there maybe multiple 
    columns contains salmon quant files. This notebook identifies the correct files
 * sample output for TCGA_READ
   ```
   $ pwd
   /Users/andrewdavidson/googleUCSC/kimLab/terraDataModels
   (base) $ ls test-aedavids-proj/TCGA/TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab/generateTCGAMatrixCreationScripts.ipynb.out
   TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_quantFile.csv
   TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_105_copyFromTerraToLocal.sh
   TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_105_copyFromTerraToNativeGCP.sh*
   TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_colData.csv
   
   $ head -n 2 TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab_105_copyFromTerraToLocal.sh
   mkdir -p ./TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab/
   gsutil -m cp gs://fc-secure-8a69fc00-b6c9-4179-aee5-f1e47a4475dd/34b2bbfb-4f9a-41d4-bfd8-b55a8e1987de/quantify/ef52a514-2fc0-4e85-946a-2bbbbc56ab96/call-salmon_paired_reads/READ-AF-2687-TP.quant.sf.gz ./TCGA_READ_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab/
   
   ```

# TODO:
- 1) ~~clean up kimLab/terraDataModels git repo~~
- 2) ~~backup / download TCGA data models~~
  * see kimLab/terraDataModels/R/downloadTCGA.Rmd
  * line 126 failed workspace='TCGA_DLBC_ControlledAccess_V1-0_DATA _edu_ucsc_kim_lab'
    + https://rawls.dsde-prod.broadinstitute.org/api/workspaces/test-aedavids-proj/TCGA_DLBC_ControlledAccess_V1-0_DATA _edu_ucsc_kim_lab/entities/participant
    ```
    ERROR: Error: 'avtable' failed:
    Internal Server Error (HTTP 500).
    Illegal URI reference: Invalid input ' ', expected '/', 'EOI', '#', '?' or pchar (line 1, column 114): https://rawls.dsde-prod.broadinstitute.org/api/workspaces/test-aedavids-proj/TCGA_DLBC_ControlledAccess_V1-0_DATA _edu_ucsc_kim_lab/entities/participant
   ```
   
- 3) review findSamplesWithQuantFile() in terra generateTCGAMatrixCreationScripts.ipynb
- 4) create notebook that generates delete script
  * for each tcga project
    + list cols with 'aux' or 'quant'
    + select cols that do not '3' in name
    write script
- 5 test delete script on smaller TCGA workspace
- 6) delete bad cols in data model
- 7) clean all
- 8) delete bad cols in data models
- 9) create backup using  kimLab/terraDataModels/R/downloadTCGA.Rmd
- 10) delete persistent disks
  * see [terra clusters](https://app.terra.bio/?utm_source=Terra+App&utm_campaign=a94b7be125-PD_PPW_migration_Campaign&utm_medium=email&utm_term=0_ea2ec28eda-a94b7be125-1366048874&ct=t(PD_PPW_migration_campaign)&mc_cid=a94b7be125&mc_eid=9a6914636f#clusters)
- 11) do we want to rewrite the code that combins the salmon quant files into a single matrix.tsv? Use wdl scatterAndGather?
