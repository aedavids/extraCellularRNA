# How to download data from exRNA.org


## Tips for managining bulk downloads
- [https://webhostinghero.org/how-to-create-a-process-group-in-linux/](https://webhostinghero.org/how-to-create-a-process-group-in-linux/)
- [https://www.howtoforge.com/linux-pstree-command/](https://www.howtoforge.com/linux-pstree-command/)
- use setsid 
This will let your batch job continue running after you log out.
setsid will make it easy to kill your bulk download process and any child process. 
See batchDownloadExRNA.orgData.sh for example. There are some trick for passing arguments 
to the child processes
```
$ setsid sh -c 'set -x; myBatchJob.sh ' > $scriptLog 2>&1 &
```

- monitoring your batch jobs using log files
if redirected stdout and stderr as in the above you can use the tail -f command

- us2 ps to montor your jobs
```
[aedavids@plaza bin]$ ps -e -o pid,ppid,pgid,command,user |head -n 1; ps -e -o pid,ppid,pgid,command,user |grep aedavids
 PID  PPID  PGID COMMAND                     USER
22866  3386 22866 sshd: aedavids [priv]       root
22867  3386 22867 sshd: aedavids [priv]       root
22871 22866 22866 sshd: aedavids@pts/84       aedavids
22872 22867 22867 sshd: aedavids@pts/85       aedavids
22873 22871 22873 -bash                       aedavids
22874 22872 22874 -bash                       aedavids
26644     1 26644 sh -c set -x; parseInstruct aedavids
26645 26644 26644 /bin/bash ./parseInstructio aedavids
26727 26645 26644 /bin/bash ./downLoadExRNA.o aedavids
26730 26727 26644 sleep 10                    aedavids
26739 22874 26739 ps -e -o pid,ppid,pgid,comm aedavids
26740 22874 26739 grep --color=auto aedavids  aedavids
```

- use pstree
```
(base) [aedavids@plaza bin]$ pstree aedavids
sshd───bash

sshd───bash───pstree
```
- killing all the childred
```
kill -TERM -[pgid]
```

## Down load instructions
- FASTQ files are big. You may want to test using the meta files. They are much smaller
- you can always re run the downloads. The scripts will skip files that have been previously down loaded. 
```
$ ls -1 extraCellularRNA/bin/exRNADownload/test*
testBioSampleDownLoadInsr.aa
testDonorMetaDownLoadInsr.aa
testExprMetaDownLoadInsr.aa
testFASTQDownLoadInstr.aa
```



Check the log files you should a "warning download skipped"

1. select a data set and follow the directions on [bulk data download tutorial](http://genboree.org/theCommons/projects/exrna-mads/wiki/Downloading_Data_and_Metadata_from_the_exRNA_Atlas)
to download the data set instructions. You will get file like extraCellularRNA/bin/exRNADownload/exRNA_FASTQ_test.tsv

2. split the instruction into smaller batch files. In the example bellow
    * sed is used to remove lines that start with #
    * cut is used to select the following columns
        + download file name
        + download url
        + sample id
        + condition
        + bio sample assession
    * split is used to break the set of instructions into a set of 10 line files. 
      These file will be in the current and be of the form ourFileprefiex.aa, ourFileprefiex.ab, ...
```
sed '/#/d' ${exRNADownLoad_tsv} | cut -f 1,2,3,4,16 | split -l 10 - ourFileprefiex.

(base) $ sed '/#/d' exRNA_FASTQ_test.tsv  | cut -f 1,2,3,4,16 | split -l 10 - testDownLoadInstr.
(base) $ ls testDownLoadInstr.aa 
testDownLoadInstr.aa

(base) $ cat testDownLoadInstr.aa 
Sample_2S14.fastq.zip   ftp://ftp.genboree.org/exRNA-atlas/grp/Extracellular%20RNA%20Atlas/db/exRNA%20Repository%20-%20hg19/file/exRNA-atlas/exceRptPipeline_v4.6.2/TPATE1-human-plasma-healthyVsCancer-2016-10-17/sample_Sample_2S14_fastq/rawInput/Sample_2S14.fastq.zip  Sample_2S14 Colon Carcinoma EXR-TPATE1NEBLIB2S14-BS
Sample_1S25.fastq.zip   ftp://ftp.genboree.org/exRNA-atlas/grp/Extracellular%20RNA%20Atlas/db/exRNA%20Repository%20-%20hg19/file/exRNA-atlas/exceRptPipeline_v4.6.2/TPATE1-human-plasma-healthyVsCancer-2016-10-17/sample_Sample_1S25_fastq/rawInput/Sample_1S25.fastq.zip  Sample_1S25 Colon Carcinoma EXR-TPATE1NEBLIB1S25-BS
Sample_N20.fastq.zip    ftp://ftp.genboree.org/exRNA-atlas/grp/Extracellular%20RNA%20Atlas/db/exRNA%20Repository%20-%20hg19/file/exRNA-atlas/exceRptPipeline_v4.6.2/TPATE1-human-plasma-healthyVsCancer-2016-10-17/sample_Sample_N20_fastq/rawInput/Sample_N20.fastq.zip    Sample_N20  Healthy Control EXR-TPATE1NEBLIBN20-BS
Sample_4S2.fastq.zip    ftp://ftp.genboree.org/exRNA-atlas/grp/Extracellular%20RNA%20Atlas/db/exRNA%20Repository%20-%20hg19/file/exRNA-atlas/exceRptPipeline_v4.6.2/TPATE1-human-plasma-healthyVsCancer-2016-10-17/sample_Sample_4S2_fastq/rawInput/Sample_4S2.fastq.zip    Sample_4S2  Colon Carcinoma EXR-TPATE1NEBLIB4S2-BS
```

3. configure mainDownload.sh
mainDownLoad.sh controls how many batches should run concurrently and what type of patch files should be run. 
When mainDownload.sh executes, it will finish quickly. All it does is start batches so that they will
continue to run after you log out. You can find the output from each batch in ./logs


