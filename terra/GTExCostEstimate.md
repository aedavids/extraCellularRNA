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

Our 'pipeline' is run as two separte workflows, SamToFastq and SalmonPairedReadQuantTask. 

## SamToFastq_copy v.1
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
  

## SamToFastq v 1 premetive should be have been cheaper
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

## SalmonPairedReadQuantTask V 17

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


# Workflow Cost Estimator Notebook

###  SamToFastq_copy v.1

submission_id = "1a0771a6-983d-47ea-8d4b-a3c2bee4aaa9"
```
                          workflow_id                      task_name  cpus  memory    disk  duration    cost
 0c4e8d80-44b5-4bb3-be51-546fa1b3eae7                     samToFastq     2    12GB    40GB     0.34h   $0.04
 00980973-94f4-4039-842f-9f414878c45c                     samToFastq     2    12GB    40GB     0.49h   $0.06
 fdf1aed9-47d0-4315-89a4-391a8ec64a9f                     samToFastq     2    12GB    40GB     1.13h   $0.14
 345079d7-2ef7-4858-9051-5bc967ada810                     samToFastq     2    12GB    40GB     1.12h   $0.14
 e05a9ae6-dd3c-4f2e-b5b9-2cc6a6882994                     samToFastq     2    12GB    40GB     0.91h   $0.11
 9a984b71-f263-4a52-9c8d-c7627363f637                     samToFastq     2    12GB    40GB     0.43h   $0.05
 2ad8f96d-69a1-4c4f-accc-9153f45fde43                     samToFastq     2    12GB    40GB     0.69h   $0.08
 2f054303-3273-4372-a43f-edc81e0932f8                     samToFastq     2    12GB    40GB     0.68h   $0.08
 08a69ff6-6a8c-4fab-a729-3ccf3a9b7d83                     samToFastq     2    12GB    40GB     0.41h   $0.05
 378a7beb-5e79-48c2-a97a-87c3c111a79e                     samToFastq     2    12GB    40GB     0.52h   $0.06
                                                                                           total_cost: $0.82
```

### SamToFastq v 1 preembtiple containers should be have been cheaper

Submission a5f18783-cb3f-4db5-8b9f-ff68cc0b59d2
```
                          workflow_id                      task_name  cpus  memory    disk  duration    cost
 42051ce2-cc30-47e1-a328-c36ddc34201e                     samToFastq     2    12GB    40GB     0.04h   $0.00
 42051ce2-cc30-47e1-a328-c36ddc34201e                     samToFastq     2    12GB    40GB     1.01h   $0.03
 7f8516d4-f5e2-4dde-9459-165e0a7af3a4                     samToFastq     2    12GB    40GB     1.21h   $0.03
 07834be0-c24e-4c4e-815c-3d712a57da10                     samToFastq     2    12GB    40GB     0.63h   $0.02
 d8e8ad2f-fb25-4c63-bced-c056652ba3bf                     samToFastq     2    12GB    40GB     0.99h   $0.03
 ab01af17-47d4-4e2c-9ed6-206bd668451c                     samToFastq     2    12GB    40GB     0.61h   $0.02
 d86185fa-9358-4c5a-9814-85a48f30688a                     samToFastq     2    12GB    40GB     0.64h   $0.02
 c1fc415a-92c7-48bb-b90b-c8e83fba043b                     samToFastq     2    12GB    40GB     0.41h   $0.01
 1ffb41ca-c745-4e54-bd7d-539a524d321a                     samToFastq     2    12GB    40GB     0.37h   $0.01
 b8096bc8-bc8f-41b3-b20f-35274c52bda5                     samToFastq     2    12GB    40GB     0.67h   $0.02
 976c6b30-4417-4902-8f58-19cf4f2663e1                     samToFastq     2    12GB    40GB     0.54h   $0.01
                                                                                           total_cost: $0.20
```

### SalmonPairedReadQuantTask V 17

Submission a3746d28-710d-4d10-b1d2-dc34dfa3c6fc
```
                          workflow_id                      task_name  cpus  memory    disk  duration    cost
 423d0833-8c01-4e89-bc9e-05df6371180f            salmon_paired_reads    10    64GB    80GB     0.63h   $0.09
 56cd92ce-a524-408e-bf10-f68b4708709d            salmon_paired_reads    10    64GB    80GB     1.86h   $0.25
 0318704e-8adb-40ac-9064-a5dc1c404c50            salmon_paired_reads    10    64GB    80GB     2.01h   $0.27
 0318704e-8adb-40ac-9064-a5dc1c404c50            salmon_paired_reads    10    64GB    80GB     0.76h   $0.10
 55e98c87-baca-455e-87be-b141b159c5c8            salmon_paired_reads    10    64GB    80GB     1.12h   $0.15
 64884339-abfe-41d8-84e7-5dac970136bc            salmon_paired_reads    10    64GB    80GB     2.70h   $0.36
 8220a3ad-42d3-463c-925a-7e780e6a5269            salmon_paired_reads    10    64GB    80GB     3.39h   $0.46
 7f8ca915-79db-4f63-bd5e-44519d79ce12            salmon_paired_reads    10    64GB    80GB     4.08h   $0.55
 d3da22b7-9ca0-4e79-a9b2-ce79d71a5ad4            salmon_paired_reads    10    64GB    80GB     3.46h   $0.47
 c8b423fc-9c1c-4a3e-9a9d-1562fc29a60f            salmon_paired_reads    10    64GB    80GB     0.53h   $0.07
 9ee645ee-a247-4fca-be83-ef357bd498d5            salmon_paired_reads    10    64GB    80GB     1.37h   $0.18
                                                                                           total_cost: $2.94
```


### Explore costs for potential workflow configurations and runtimes.
I think this is the cost for a single docker container

```
 cpus   memory     disk  runtime  preemptible     cost
      10     64GB    700GB       5h        False    $3.27
       8     32GB    700GB      10h        False    $4.46
      10     64GB    700GB       5h         True    $0.84
       8     32GB    700GB      10h         True    $1.24
       8     32GB    400GB      10h         True    $1.08
       8     32GB    100GB      10h         True    $0.91
```
