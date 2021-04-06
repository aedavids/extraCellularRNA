# GTEx app.terra.bio cost estimate

All test where run in a clone of the the AnVIL workspace
[https://app.terra.bio/#workspaces/test-aedavids-proj/AnVIL_GTEx_V8_hg38_AEDWIP](https://app.terra.bio/#workspaces/test-aedavids-proj/AnVIL_GTEx_V8_hg38_AEDWIP)



I select 10 of hte 328 Pancreas samples. Terra runs each of these in a separate docker container

```
GTEX-11GSP-0426-SM-5A5KX (sample)
GTEX-11I78-0626-SM-5A5LZ (sample)
GTEX-11LCK-0226-SM-5A5M6 (sample)
GTEX-11NSD-0526-SM-5A5LT (sample)
GTEX-11ONC-0526-SM-5BC57 (sample)
GTEX-11P7K-0526-SM-5BC5I (sample)
GTEX-11TT1-0326-SM-5LUAY (sample)
GTEX-11VI4-0426-SM-5EGHZ (sample)
GTEX-11XUK-0626-SM-5N9ES (sample)
GTEX-1211K-1126-SM-5EGGB (sample)
```

Our 'pipeline' is run as two separte workflows. 

1. SamToFastq_copy v.1
   - Submission 1a0771a6-983d-47ea-8d4b-a3c2bee4aaa9
   - panc-10Samples-2021-04-06T00-57-33 sample_set
   - start Apr 6, 2021, 8:25 AM
   - 10 gb of memory was not enough to run all the samples
   ```
   {
   "samToFastqTest.samToFastq.diskSpaceGb": "${40}",
   "samToFastqTest.samToFastq.inputBam": "${this.bam_file}",
   "samToFastqTest.samToFastq.memoryGb": "${12}",
   "samToFastqTest.samToFastq.sampleName": "${this.sample_id}"
   }
   ```


premetive would have been cheaper


storage costs?

egress costs?

