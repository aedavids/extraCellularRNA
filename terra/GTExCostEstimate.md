# GTEx app.terra.bio cost estimate
April 2021 
Andy Davidson
aedavids@ucsc.edu

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
   - end   Apr 6, 2021, 9:34 AM
   - last batch completed on 
   - 10 gb of memory was not enough to run all the samples
   ```
   {
   "samToFastqTest.samToFastq.diskSpaceGb": "${40}",
   "samToFastqTest.samToFastq.inputBam": "${this.bam_file}",
   "samToFastqTest.samToFastq.memoryGb": "${12}",
   "samToFastqTest.samToFastq.sampleName": "${this.sample_id}"
   }
   ```
   
   output wound up in a weird output column
   ```
   this. unmapped_first_end_fastq	
   this.unmapped_second_end_fastq	
   this.unmapped_unpairedFastq
   ```
  

2. SamToFastq v 1 premetive should be have been cheaper
   - Submission a5f18783-cb3f-4db5-8b9f-ff68cc0b59d2
   - panc-10Samples-2021-04-06T00-57-33
   - Source: aedavids.ucsc.edu/samToFastq/1
   - preemptible is set to 3 by default
   - start Apr 6, 2021, 10:40 AM
   - end   Apr 6, 2021, 11:53 AM
   - output
   ```
   this.firstEndFastq
   this.secondEndFastq
   this.unpairedFastq
   ```


long term storage costs?

egress costs? Worflow says 'Free egress to this data is available' Not clear if we are being charge and if so what we have to do to get free access

3. SalmonPairedReadQuantTask V 17

    - Submission a3746d28-710d-4d10-b1d2-dc34dfa3c6fc
    - workflow config https://app.terra.bio/#workspaces/test-aedavids-proj/AnVIL_GTEx_V8_hg38_AEDWIP/workflows/aedavids.ucsc.edu/SalmonPairedReadQuantTask
    - start Apr 6, 2021, 12:06 PM
    - end   Apr 6, 2021, 4:13 PM
    - ref: sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx.tar.gz
    - used default parameters in workflow. required a suprising amount of memory
    ```
    String dockerImg = 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
    Int runtime_cpu = 8
    Int memoryGb = 64
    Int diskSpaceGb = 80
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    Int runtime_preemptible = 3
    ```

