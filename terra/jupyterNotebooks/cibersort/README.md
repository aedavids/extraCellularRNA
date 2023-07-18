# terra/jupyterNotebooks/cibersort

notebooks are listed in chronologic order

**overview**

- CIBERSORTx-Tutorial.ipynb
  * create a trivial sample data set that demonstrates 
      + how to format gene signature and mixture matrix files
      + example of cibersort output
      
 - CreateGeneSignatureMatrixOverview.ipynb
   * illustrates how we construct the gene signature matrix using a small data set
   
- createCiberSortGeneSignatureMatrix.ipynb
  * given results files extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb and results form 1vsAll analysis, creates  a gene signature matrix
  
- createFractionMatrixOverview.ipynb
  * explians how we create the expectedFractions matrix in createCibersortMixtureMatrix

- createCibersortMixtureMatrix.ipynb
    * each sample/row in the mixture model is a single type
    * expected fractions matrix
      - can be used to create other test mixture models
        * by gender
        * by data set type (GTEx or TCGA)
        * random combinations
        * ...
    * random mixture matrix
      - all info removed, can be used to create a base line to compare models agains
      

- closedFormDeconvolutionTest.ipynb
  * cibersort docker took 3.5 days to run GTEx_TCGA training data set. This notebook explores using a closed form solution. It does not appear to work. There are some comments about how to improve it.


- preliminaryCibersortResults.ipynb
  * explore results from inital run on GTEx_TCGA training data set using extraCellularRNA/terra/deseq/data/private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design:~gender+category-padj:0.001-lfc:2.0-n:25/ciberSort
  * plots, wisker_box, scatter, ... fractions table, ...

- evaluateParallelTestResults.ipynb
  * cibersort took 3.5 days to run GTEx_TCGA training data set. This notebook attempts to reduce compute time by spliting mixture into pars and running dockers in parrallel
    * <span style="color:green">It is possible to split traning data set in to shards run cibersort concurrently and re-assemble results</span>
    * ref: extraCellularRNA/terra/cibersortx/cibersortParallelization.md

- clusters
  * preliminaryCibersortClusterPlots.ipynb TODO louvain + UMAP
  * preliminaryResultsCibersortFractionsPCA.ipynb
    - PCA plot male/female, and GTEx/TCGA
  * preliminaryResultsCibersortFractionsUMAP.ipynb
    - uses shapes and colors to show how all 83 classes cluster

- volcanoPlots.ipynb
  * used a notebook to keep track of how all the plots where created
  
- contingencyTables.ipynb
  * also known as a cross tabulation or crosstab
  * can we use contigency tables to visuallize deconvolution data
  * the fractional tables in  preliminaryCibersortResults.ipynb are hard to see
  * uses indicator variable to show kidney samples might have traces of cancer signals
  * pivottables is implement using react.js . You can click on 'pop out' to pop out a browser window with an html page you can cut and pages the img from. The argument outfile_path does not work the way you might like, it create an html page and saves it the computer the browser is running on. You can not save the server
  
- fractionsAsMulticlassClassificationPOC.ipynb
  * sample code used to figure out how to do the anlysis
  * contains lots of debug 
  
- fractionsAsMulticlassClassification.ipynb
  * analysize as if fractions where output layer of multclass classifier
  * this is a copy of fractionsAsMulticlassClassificationPOC.ipynb with out all the debug
  * output
    - a series of confusion matrixes. There are two main flavors. 
      + actual TP and FP counts
      + row percentages, it of all the sample of type "A" what percent where TP or call as type B, C, ..
      
  
  
- verifyCIBERSORTxFractionsWorkflow.ipynb 3/28/23
  * verify scatter/gather produces same results as run with entire mixture file

- 7/18/2023 scalingBug/README.md
  * I wanted to save copies of these notebooks. I did not want to merge these notebooks back to main are
  