# keep track of where our containers came from

- rocker/rstudio:4.0.0-ubuntu18.04
    + aedavids/ggplot2
    
- rocker/rstudio:3.5.0 Debian GNU/Linux 9
    + aedavids/biocworkshops
    + <span style="color:red">fail can not install DeSeq2</span>
   
- bioconductor/bioconductor_docker:devel Ubuntu 18.04.4 LTS R4.0.0
    + aedavids/biocworkshop2018desq2
        * installed packages for desq2 workshop using BiocManager::install()y
    +  aedavids/extra_cellular_rna_broken
       * installed tidyr in rStudio server console

-  bioconductor/bioconductor_docker:devel "BIOCONDUCTOR_DOCKER_VERSION=3.12.31"
   <span style="color:red">I do not think this image is portable. lots of files get compiled</>
    run following R commands in RStudio console
    ```
    if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")

    BiocManager::install("DESeq2")
    BiocManager::install("apeglm")
    
    install.packages("tidyr")
    install.packages("readr")
    install.packages("dplyr")
    install.packages("rjson")
    install.packages("rmarkdown")
    install.packages("tximport")
    ```
