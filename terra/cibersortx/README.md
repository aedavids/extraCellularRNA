# How to use Cibersort

```
Andrew E. Davidson
aedavids@ucsc.edu
```

ref:
- [https://cibersortx.stanford.edu/download.php](https://cibersortx.stanford.edu/download.php) select Access instructions (zip)


## cibersort input data

long term storage: /private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/{best,up}/*/ciberSort/

created by 
- extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
- extraCellularRNA/terra/jupyterNotebooks/cibersort/createCibersortMixtureMatrix.ipynb 

## Running fractions analysis

### Test 1: use small data set to debug docker arguments
use small hand craft test data. see extraCellularRNA/terra/jupyterNotebooks/cibersort/CIBERSORTx-Tutorial.ipynb . In this notebook we uploaded the hand crafted data set and ran on stanford's server. We should get the same resuls when we run using docker localling. Copy data generated input files from mac
```
cd /Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/terra/jupyterNotebooks/cibersort
ssh mustard mkdir /scratch/aedavids/GTEx_TCGA/geneSignatureProfiles/test

scp mixture.txt signatureGenes.txt mustard:/scratch/aedavids/GTEx_TCGA/geneSignatureProfiles/test/
mixture.txt                                         100%  246    20.8KB/s   00:00    
signatureGenes.txt                                  100%  134    11.2KB/s   00:00   

cat signatureGenes.txt
name	T1	T2	T3
G1	1.0	0.0	0.0
G2	1.0	0.0	0.0
G3	1.0	1.0	0.0
G4	1.0	1.0	0.0
G5	0.0	1.0	1.0
G6	0.0	0.0	1.0
G7	0.0	0.0	0.0
G8	0.0	1.0	0.0

cat mixture.txt
sampleTitle	S1	S2	S3	S4	S5	S6
G1	1.0	0.0	0.0	1.0	1.0	0.0
G2	1.0	0.0	0.0	1.0	1.0	0.0
G3	1.0	1.0	0.0	2.0	1.0	1.0
G4	1.0	1.0	0.0	2.0	1.0	1.0
G5	0.0	1.0	1.0	1.0	1.0	2.0
G6	0.0	0.0	1.0	0.0	1.0	1.0
G7	0.0	0.0	0.0	0.0	0.0	0.0
G8	0.0	1.0	0.0	1.0	0.0	1.0
```

configure and run
```
dst=/scratch/aedavids/GTEx_TCGA/geneSignatureProfiles/test
cd ${dst}
mixtureMatrix=${dst}/mixture.txt
signatureMatrix=${dst}/signatureGenes.txt

# mount points must be full paths
docker run \
    -v `pwd`:/src/data \
    -v ${outDir}:/src/outdir cibersortx/fractions \
    --username aedavids@ucsc.edu \
    --token 3f561ab6d4cf373d11f23d8e205b4b72 \
    --mixture ${mixtureMatrix}\
    --sigmatrix ${signatureMatrix}\
    --perm 100 \
    --label $jobId \
    --QN FALSE \
    --verbose TRUE
```


### Test 2: run docker with real data

extraCellularRNA/terra/cibersortx/bin/run_cibersortx_fractions.sh. This script will run docker detached and will automatically remove container on exit. No need to use setsid

output files

```
cd /scratch/aedavids/cibersort.out
./ciberSort.debug.out.2022-10-18-07.37.32-PDT 
./2022-09-27-simpleTest/CIBERSORTx_2022-09-27-simpleTest_Results.txt
./GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT/CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt
```

### results from run Crete (Chania)  sat 10/22/22

```
$ ll /scratch/aedavids/cibersort.out/GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT/CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt 
-rw-r--r-- 1 root root 27M Oct 21 18:42 /scratch/aedavids/cibersort.out/GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT/CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt
```

input was training data set. each sample was from a single type. ie lung, PAAD, ...


approximate total run time 3.5 days
```
  2022-10-18-07.40.54-PDT
-      Oct 21 18:42
-----------------------

 00:20 # 10/18 08:00
+16:00 # 10/19
+24:00 # 10/20
+24:00 # 10/21
+18:42
---------
 83:02

```

### explore results file
15801 samples
```
$ wc -l *Results.txt
15802 CIBERSORTx_GTEx_TCGA_TrainGroupby_mixture-2022-10-18-07.40.54-PDT_Results.txt
```

number of columns in results
```
$ ls *Results.txt > t
$ sft
$ head -n 1 $f | sed -e 's/\t/\n/g' |wc -l
87
```

These samples are from GTEx. looking at the *Results.txt file alone, it is hard to know what the true tissue type is
Correlation looks good. we expect that samples only deconvole into a single type: TODO write code to check if fraction match expection

```
$ head  $f | cut -f 1,2,3,84,85,86,87
Mixture	ACC	Adipose_Subcutaneous	Whole_Blood	P-value	Correlation	RMSE
GTEX-1117F-0226-SM-5GZZ7	0.00000238818947786	0.32088468153413718	0.00000000000000000	0.00000000000000000	0.98542467644677711	0.92645530582020463
GTEX-1117F-0526-SM-5EGHJ	0.00002175019341662	0.29529977326057638	0.00000000000000000	0.00000000000000000	0.97979080799155893	0.93469518567034882
GTEX-1117F-0726-SM-5GIEN	0.00015259963369839	0.05445776007930397	0.00266438319091741	0.00000000000000000	0.98490561378709340	0.44791559959147931
GTEX-1117F-2826-SM-5GZXL	0.00000000000000000	0.29004114942525649	0.01183538119295482	0.00000000000000000	0.98846407133779690	0.90767508784039519
GTEX-1117F-3226-SM-5N9CT	0.00021702689123660	0.00445898085087070	0.00043860247181520	0.00000000000000000	0.91705227835834124	0.94955257328411946
GTEX-111CU-0326-SM-5GZXO	0.00012839158257415	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.96999619583782004	0.71740096723763391
GTEX-111CU-0426-SM-5GZY1	0.00017911589602714	0.06144797178936864	0.00000000000000000	0.00000000000000000	0.96897428384635731	0.80545057475847193
GTEX-111CU-0526-SM-5EGHK	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.90218261892042595	1.23410866579020806
GTEX-111CU-0626-SM-5EGHL	0.00073656971878856	0.20411159050680064	0.00000000000000000	0.00000000000000000	0.82951499979038523	0.95537827265096409
```

these samples are from TCGA. Looking at the results.txt file alone, it is hard to know what the true cancer type is 
correlation for some of the samples does not look great. Keep in mind we have not done any tunning yet.
```
$ tail  $f | cut -f 1,2,3,84,85,86,87
UVM-WC-A881-TP	0.00000863062071495	0.00309889939958328	0.00000000000000000	0.00000000000000000	0.79120013970914405	0.95817331960717933
UVM-WC-A882-TP	0.00000000000000000	0.00000000000000000	0.01065985495150996	0.00000000000000000	0.64824208271638017	0.95821000145648949
UVM-WC-A883-TP	0.00002268063872719	0.00196036493385387	0.08205200129025236	0.00000000000000000	0.93931661338215500	0.71381242623382180
UVM-WC-A884-TP	0.00015138225222482	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.76022553352495725	0.96080226819233827
UVM-WC-A888-TP	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.82092837581099531	0.96390208810991762
UVM-WC-AA9E-TP	0.00000000000000000	0.00000000000000000	0.00677352970766185	0.00000000000000000	0.70668108537098062	0.95809575383809553
UVM-YZ-A980-TP	0.00095955758606998	0.00314759591266955	0.00000000000000000	0.00000000000000000	0.63927936917453354	0.97024860638613908
UVM-YZ-A983-TP	0.00000894227722533	0.00000000000000000	0.00000000000000000	0.01000000000000001	0.61636688808242612	0.96229386850789200
UVM-YZ-A984-TP	0.00000000000000000	0.00000000000000000	0.00312798052610109	0.00000000000000000	0.79698294355561805	0.95684729590763562
UVM-YZ-A985-TP	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.00000000000000000	0.91638874519350044	0.95127364087000732
```

Our untuned model can only account for 38% of the variation in sample UVM-YZ-A983-TP. UVM is Uveal Melanoma

```
>>> rSquared = 0.616 ** 2
>>> rSquared
0.379456
```
