# ref https://api.firecloud.org/ga4gh/v1/tools/aedavids:SamToFastq_copy/versions/1/plain-WDL/descriptor
# Andrew Davidson
# aedavids@ucsc.edu
workflow samToFastqTest {
    call samToFastq
}

task samToFastq {
    File inputBam
    String sampleName
    Int memoryGb
    Int diskSpaceGb

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

    command <<<
        # use cpuinfo to debug preemptiple/cpu/threading quotas and performance
        cat /proc/cpuinfo
        
	java -jar /usr/gitc/picard.jar SamToFastq \
	TMP_DIR=. \
	INPUT=${inputBam} \
	FASTQ=${sampleName}.1.fastq.gz \
	INTERLEAVE=false \
	SECOND_END_FASTQ=${sampleName}.2.fastq.gz \
	INCLUDE_NON_PF_READS=true \
	CLIPPING_ATTRIBUTE=XT \
	CLIPPING_ACTION=2 \
	UNPAIRED_FASTQ=${sampleName}.unpaired.fastq.gz
    >>>

    output {
	File firstEndFastq = "${sampleName}.1.fastq.gz"
	File secondEndFastq = "${sampleName}.2.fastq.gz"
	File unpairedFastq = "${sampleName}.unpaired.fastq.gz"
    }

    runtime {
	docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
	memory: "${memoryGb} GB"
	cpu: "1"
	disks: "local-disk ${diskSpaceGb} HDD"
        preemptible: '${runtime_preemptible}' 
    }
}
