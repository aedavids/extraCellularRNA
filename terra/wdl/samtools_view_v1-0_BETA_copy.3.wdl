task samtools_view {

    File bam_file
    #File bam_index
    #String prefix
    String sampleName
    String? options # -b -f 4 
    #String? region
    ##File? reference_fasta
    # File? reference_fasta_index

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running 'samtools view'.")
        samtools view ${options} ${bam_file}  > ${sampleName}.umap.bam
        echo $(date +"[%b %d %H:%M:%S] done.")
    }

    output {
        File unmappedBAM = "${sampleName}.umap.bam"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow samtools_view_workflow {
    call samtools_view
}
