# Test terra caching
## Abstract
<span style="color:red">almost works as expected.</span> 

The second run ran in 5 mins. 
No new work was done however the results from the previous run where copied into the cache of the second. 
The sample table has the cached bucket URLS from the second run. The original bam input file bucket
URLS remain unchanged

## results
cromwell searched the cache of previous run jobs for ones with exact same command and inputs. If found will use results from previous run.

ref: Notebook p. 167, 5/4/2021 "Terra Robust Workflows"

## Test step 1 results JobHistory

```
ID:
5290da5f-5fe3-4e66-b3a1-9262d2e1ad87
workspace-id: f5aa8a37-78e5-45f6-9c59-c643016f7d97submission-id: eb866962-a367-4463-a1b2-77c069150749

duration: 1 hr 7 min:

sample 1
"links" -> job manager
https://job-manager.dsde-prod.broadinstitute.org/jobs/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87

firstEndFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87/call-bamToFastq/GTEX-111CU-0526-SM-5EGHK.1.fastq.gz
secondEndFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87/call-bamToFastq/GTEX-111CU-0526-SM-5EGHK.2.fastq.gz
unpairedFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87/call-bamToFastq/GTEX-111CU-0526-SM-5EGHK.unpaired.fastq.gz


aux_info
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87/call-salmon_paired_reads/GTEX-111CU-0526-SM-5EGHK.aux_info.tar.gz
quantFile
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87/call-salmon_paired_reads/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz

workflow dashboard
quantify.bamToFastq call chaching result : cache Miss 
quantify.salmon_paired_reads: call chaching result : cache Miss 
https://app.terra.bio/#workspaces/test-aedavids-proj/AnVIL_GTEx_V8_hg38_edu_ucsc_kim_lab/job_history/eb866962-a367-4463-a1b2-77c069150749/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87

sample 2
https://job-manager.dsde-prod.broadinstitute.org/jobs/6696a449-2d36-4a1b-96ba-91c443eb93a1


firstEndFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/6696a449-2d36-4a1b-96ba-91c443eb93a1/call-bamToFastq/GTEX-111YS-1226-SM-5EGGJ.1.fastq.gz
secondEndFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/6696a449-2d36-4a1b-96ba-91c443eb93a1/call-bamToFastq/GTEX-111YS-1226-SM-5EGGJ.2.fastq.gz
unpairedFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/6696a449-2d36-4a1b-96ba-91c443eb93a1/call-bamToFastq/GTEX-111YS-1226-SM-5EGGJ.unpaired.fastq.gz


aux_info
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/6696a449-2d36-4a1b-96ba-91c443eb93a1/call-salmon_paired_reads/GTEX-111YS-1226-SM-5EGGJ.aux_info.tar.gz
quantFile
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/6696a449-2d36-4a1b-96ba-91c443eb93a1/call-salmon_paired_reads/GTEX-111YS-1226-SM-5EGGJ.quant.sf.gz


## Test step2
rerun. DO NOT CHANGE INPUT OR OUTPUT

expected result is nothing was recomputed.

submission ID 3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b
duration: 5 min

sample 1 list output
https://job-manager.dsde-prod.broadinstitute.org/jobs/e521dbca-0912-4269-9d22-4f12a78dbd25

firstEndFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-bamToFastq/cacheCopy/GTEX-111CU-0526-SM-5EGHK.1.fastq.gz
secondEndFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-bamToFastq/cacheCopy/GTEX-111CU-0526-SM-5EGHK.2.fastq.gz
unpairedFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-bamToFastq/cacheCopy/GTEX-111CU-0526-SM-5EGHK.unpaired.fastq.gz


aux_info
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-salmon_paired_reads/cacheCopy/GTEX-111CU-0526-SM-5EGHK.aux_info.tar.gz
quantFile
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-salmon_paired_reads/cacheCopy/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz

workflow dash board
https://app.terra.bio/#workspaces/test-aedavids-proj/AnVIL_GTEx_V8_hg38_edu_ucsc_kim_lab/job_history/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/e521dbca-0912-4269-9d22-4f12a78dbd25

quantify.bamToFastq: calling Cache Result Cache Hit: 5290da5f-

sample 1 list output
https://job-manager.dsde-prod.broadinstitute.org/jobs/6b203348-2dd9-443c-847f-6be762993b74


firstEndFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/6b203348-2dd9-443c-847f-6be762993b74/call-bamToFastq/cacheCopy/GTEX-111YS-1226-SM-5EGGJ.1.fastq.gz
secondEndFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/6b203348-2dd9-443c-847f-6be762993b74/call-bamToFastq/cacheCopy/GTEX-111YS-1226-SM-5EGGJ.2.fastq.gz
unpairedFastq
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/6b203348-2dd9-443c-847f-6be762993b74/call-bamToFastq/cacheCopy/GTEX-111YS-1226-SM-5EGGJ.unpaired.fastq.gz


aux_info
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/6b203348-2dd9-443c-847f-6be762993b74/call-salmon_paired_reads/cacheCopy/GTEX-111YS-1226-SM-5EGGJ.aux_info.tar.gz
quantFile
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/6b203348-2dd9-443c-847f-6be762993b74/call-salmon_paired_reads/cacheCopy/GTEX-111YS-1226-SM-5EGGJ.quant.sf.gz

workflow dashboard
https://app.terra.bio/#workspaces/test-aedavids-proj/AnVIL_GTEx_V8_hg38_edu_ucsc_kim_lab/job_history/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/6b203348-2dd9-443c-847f-6be762993b74

quantify.bamToFastq : Call Caching Result Cache Hit: 6696a449-2d36-4a1b-96ba-91c443eb93a1:quantify.bamToFastq:-1
quantify.salmon_paired_reads : calling caching Result: Cache Hit: 6696a449-2d36-4a1b-96ba-91c443eb93a1:quantify.bamToFastq:-1
```


## Is are the output buckets the same?
<span style="color:red">NO! it looks like the output from the first get copied into the cache of the second. We have 2 copies</span>

```
first run:
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87/call-salmon_paired_reads/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz

second run :
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-salmon_paired_reads/cacheCopy/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz
```


```
(extraCellularRNA) [aedavids@mustard extraCellularRNA]$ gsutil ls -L gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87/call-salmon_paired_reads/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz 
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/eb866962-a367-4463-a1b2-77c069150749/quantify/5290da5f-5fe3-4e66-b3a1-9262d2e1ad87/call-salmon_paired_reads/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz:
    Creation time:          Tue, 04 May 2021 18:36:40 GMT
    Update time:            Tue, 04 May 2021 18:36:40 GMT
    Storage class:          STANDARD
    Content-Length:         85468043
    Content-Type:           application/octet-stream
    Hash (crc32c):          jDLW1Q==
    Hash (md5):             IORDev+FT8cMPMKq9izplw==
    ETag:                   COvDot/VsPACEAE=
    Generation:             1620153400336875
    Metageneration:         1
    ACL:                    []
TOTAL: 1 objects, 85468043 bytes (81.51 MiB)
(extraCellularRNA) [aedavids@mustard extraCellularRNA]$ gsutil ls -L gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-salmon_paired_reads/cacheCopy/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-salmon_paired_reads/cacheCopy/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz:
    Creation time:          Tue, 04 May 2021 19:09:17 GMT
    Update time:            Tue, 04 May 2021 19:09:17 GMT
    Storage class:          STANDARD
    Content-Length:         85468043
    Content-Type:           application/octet-stream
    Hash (crc32c):          jDLW1Q==
    Hash (md5):             IORDev+FT8cMPMKq9izplw==
    ETag:                   CJ78yYTdsPACEAE=
    Generation:             1620155357625886
    Metageneration:         1
    ACL:                    []
TOTAL: 1 objects, 85468043 bytes (81.51 MiB)
```


## check sample table bucket URLS for output

```
gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-bamToFastq/cacheCopy/GTEX-111CU-0526-SM-5EGHK.1.fastq.gz

gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-salmon_paired_reads/cacheCopy/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz
```


## check  sample table bucket URL for original bam files
```
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/GTEX-111CU-0526-SM-5EGHK.Aligned.sortedByCoord.out.patched.md.bam
```
