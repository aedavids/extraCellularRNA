# extraCellularRNA/terra/cibersortxextraCellularRNA/terra/cibersortx overview

```
Andrew E. Davidson
aedavids@ucsc.edu
3/21/23
```

**Table of Contents**

- [gettingStartedWithCIBERSORTx.md](./gettingStartedWithCIBERSORTx.md)
  * ref to juypyter notes used 
    + create test data we upload to CIBERSORTx website
    + use test data to figure out how to run CIBERSORTx docker on mustard
    + ran docker on Reall Data set. (Took over three days)
- [cibersortParallelization.md](./cibersortParallelization.md)
  * documents a juypter notebook and shell script we used to prove
  we could split a mixture file into parts. Run the parts separately, combine the results and get the same resuls as if we ran original mixture data file
- [wdlTest/README.md](./wdlTest/README.md)
  * I used this directory to write a CIBERSORTx scatter/gather WDL workflow POC.
- [wdl/READM.md](./wdl/README.md)
  * describes how I developed the CIBERSORTx fractions workflow
