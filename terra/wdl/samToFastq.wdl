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
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    Int runtime_preemptible = 3

    command <<<
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
