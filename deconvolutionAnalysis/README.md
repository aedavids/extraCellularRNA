# deconvolution Analysis 
Andrew E. Davidson  
aedavids@ucsc.edu  
9/23


The deconvolution directory  contains  code to make it easier to evaluate our deconvolution results and to tune our pipline hyper parameters

deconvolutionAnalysis/bin/best25GTEx_TCGA.sh demonstrates how to run the pipeline

the pipeline does the following

1. activate the extraCellularRNA conda environment
2. run pipeline.upstreamPipeline
   a. user passes in a python module that implements a function named 
   ```
   def createSignatureGeneConfig(outDirPath : str, 
       vargs : list = None ) -> SignatureGeneConfiguration
   ```
   
   the users module implements a class derived from SignatureGeneConfiguration taht selects genes of interest. vargs are used to allow the user to pass through what ever arguments their implementation requires
   
   analysis.createBestCreateSignatureGeneConfig is a good example.
   
   b. step 1: 
       - create user defined SignatureGeneConfiguration obj
       - run SignatureGeneConfiguration.findGenes() on all candidate results file
   c. step 2:
       - create upsets plotss
   d. step 3 create a signature matrix in CIBERSORTx format
   e. step 4 create in  CIBERSORTx format
       - mixture matrix
       - expected fractions. Assume each sample is composed of single type/categor/tissue type/cancer type
       - randomized mixture matrix
3. runs CIBERSORTxFractionsWorkflow.wdl.
4. runs analysis.metrics module
   this will output standard k-way classifier stats. f-score, sensitivity, specificity, ...

## to run the unit test

run all unit test
```
conda activate extraCellularRNA
cd deconvolutionAnalysis/python/
python -m unittest discover .
```

to run a single test case 
```
# run  testFindDegrees() case in the TestUpsetPlots class
conda activate extraCellularRNA
cd deconvolutionAnalysis/python/
python plots/test/testUpsetPlots.py  TestUpsetPlots.testFindDegrees
```
