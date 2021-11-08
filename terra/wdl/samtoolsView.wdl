task samtools_view {

    File bam_file
    #File bam_index
    #String prefix
    String sampleName
    String? options # "-b -f 4"
    #String? region
    ##File? reference_fasta
    # File? reference_fasta_index

    Int memoryGb = 10
    Int diskSpaceGb = 40

    #
    # best practice 2 core min: one for os one for work
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#cpu
    # In Google Cloud: this is interpreted as "the minimum number of cores to use."
    #        
    Int runtime_cpu = 2

    #
    # https://cloud.google.com/kubernetes-engine/docs/how-to/preemptible-vms
    # instances that last a maximum of 24 hours in general, and provide no availability guarantees.
    # Preemptible VMs are priced lower than standard Compute Engine
    #
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    # Take an Int as a value that indicates the maximum number of times Cromwell should request a
    # preemptible machine for this task before defaulting back to a non-preemptible one.
    # default value: 0
    Int runtime_preemptible = 3    

    command {
        set -euo pipefail

        # use cpuinfo to debug preemptiple/cpu/threading quotas and performance
        cat /proc/cpuinfo
        
        echo $(date +"[%b %d %H:%M:%S] Running 'samtools view'.")
        samtools view ${options} ${bam_file}  > ${sampleName}.umap.bam
        echo $(date +"[%b %d %H:%M:%S] done.")
    }

    output {
        File unmappedBAM = "${sampleName}.umap.bam"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memoryGb}GB"
        disks: "local-disk ${diskSpaceGb} HDD"
        cpu: "${runtime_cpu}"
        preemptible: "${runtime_preemptible}"
    }

    meta {
        author: "Andrew Davidson, based on version by Francois Aguet"
    }
}


workflow samtools_view_workflow {
    call samtools_view
}
