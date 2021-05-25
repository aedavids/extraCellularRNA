
# note urls have version information
import "https://api.firecloud.org/ga4gh/v1/tools/aedavids.ucsc.edu:bamToFastq/versions/2/plain-WDL/descriptor" as bamToFastq
import "https://api.firecloud.org/ga4gh/v1/tools/aedavids.ucsc.edu:SalmonPairedReadQuantTask/versions/25/plain-WDL/descriptor" as salmonPairedReadQuantTask

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
    }
}
