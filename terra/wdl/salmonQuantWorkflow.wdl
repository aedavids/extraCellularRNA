
# note urls have version information
import "https://api.firecloud.org/ga4gh/v1/tools/aedavids.ucsc.edu:bamToFastq/versions/16/plain-WDL/descriptor" as bamToFastq
import "https://api.firecloud.org/ga4gh/v1/tools/aedavids.ucsc.edu:SalmonPairedReadQuantTask/versions/27/plain-WDL/descriptor" as salmonPairedReadQuantTask

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
        firstFastq=bamToFastq.firstEndFastq,
        secondFastq=bamToFastq.secondEndFastq,
    }

     output {
         # the fastq.qz files created by bamToFastq.bamToFastq are really big
         # generates a lot of storage cost. if we do not bind them to a data
         # model and we enable 'Use call caching' and
         # 'Delete intermediate outputs' then if a submission is sucesfull
         # the fastq.gz files will be deleted if the submission fail the
         # fastq.gz files will remain in the submission cache we can
         # re-run with out having to recompute
         #
         # https://support.terra.bio/hc/en-us/articles/360039681632-Saving-storage-costs-by-deleting-intermediate-files
         # https://app.terra.bio/#workspaces/help-gatk/Terra-Tools
         
         File quantFile     = salmon_paired_reads.quantFile
         File aux_info     = salmon_paired_reads.aux_info 
     }     
}
