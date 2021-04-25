
import "https://api.firecloud.org/ga4gh/v1/tools/aedavids.ucsc.edu:bamToFastq/versions/1/plain-WDL/descriptor" as bamToFastq
import "https://api.firecloud.org/ga4gh/v1/tools/aedavids.ucsc.edu:SalmonPairedReadQuantTask/versions/22/plain-WDL/descriptor" as salmonPairedReadQuantTask

workflow quantify {
    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

    String sampleId
    File inputBam
    File refIndexTarGz
    
    call bamToFastq.bamToFastq {
        input:
        inputBam=inputBam,
        sampleName=sampleId
    }

    call salmonPairedReadQuantTask.salmon_paired_reads {
        input:
        sampleId=sampleId,
        refIndexTarGz=refIndexTarGz,
        leftReads=bamToFastq.firstEndFastq,
        rightReads=bamToFastq.secondEndFastq,

        # see salmonPairedReadQuantTask.wd for more information aobut
        # choosing parameter values
        outDir="salmon.out",
        dockerImg='quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0',
        memoryGb=64,
        diskSpaceGb=80,
        runTimeCpu=8,        
        runTimePreemptible=3
    }
}
