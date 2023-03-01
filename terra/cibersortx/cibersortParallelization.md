# Can we split the mixture matrix into parts and run in parallel?
```
Andrew E. Davidson
aedavids@ucsc.edu
```

**abstract**
running cibersort on our GTEx_TCGA training data set took 3.5 days. It is possible to split training data set into shards and run cibersort concurrently

1. Create 3 mixture matric
   a. referenceMixture with 100 samples
   b. test1Mixture and test2Mixture each with 1/2 of the referenceMixture samples
2. review extraCellularRNA/terra/cibersortx/bin/run_cibersortx_fractions.sh
3. create simplified vesion
   a. three blocks. one for each mixture we want to run

# run test
- extraCellularRNA/terra/cibersortx/bin/parallelTest.sh

## evaluate results
see extraCellularRNA/terra/jupyterNotebooks/cibersort/evaluateParallelTestResults.ipynb
