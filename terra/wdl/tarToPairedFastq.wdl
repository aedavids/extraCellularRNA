#
# prepares a tar file containing fastq files so that it can be used with
# by a down stream wdl task like salmonPairedReadQuantTask.wdl
#
# the input tarfile is assume to have paired end
# fastq files. There may be replicants. The fastq files can be in subdirectories
#
# for paired fastqs two gzip files will be created. The first will contain a single file
# created by concatenating all the "1" fastq files. The second will contain all the "2" files
#
# TODO
# add support for single end fastq files
#

workflow tarToFastqTest {
    meta {
        author: "Andrew E. Davidson"
	email: "aedavids@ucsc.edu"
    }
    
    call tarToFastqTask 
}

task tarToFastqTask {

    # input
    File tarFile
    String sampleName
    Int memoryGb=32
    Int diskSpaceGb=40
    
    #
    # https://cloud.google.com/kubernetes-engine/docs/how-to/preemptible-vms
    # instances that last a maximum of 24 hours in general, and provide no availability guarantees.
    # Preemptible VMs are priced lower than standard Compute Engine
    #
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    # Take an Int as a value that indicates the maximum number of times Cromwell should request a
    # preemptible machine for this task before defaulting back to a non-preemptible one.
    # default value: 0
    Int runTimePreemptible = 3    

    command <<<
        # we do not control the execution of this script
        # cromwell/docker ? puts set -x in stderr. this makes debugging harder
        # https://ops.tips/gists/redirect-all-outputs-of-a-bash-script-to-a-file/#using-exec-to-redirect-stdout-and-stderr
        # Redirect standard error to standard out such that 
        # standard error ends up going to wherever standard
        # out goes (the file).
        #
        # this will cause the 'set -x' output to be interleaved with the output of each statement
        #
        exec 2>&1

        echo "input values"
        echo "tar: ${tarFile}"        
        echo "sampleName: ${sampleName}"

        # configure script to make debug easier
        # ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca         
        set -euxo pipefail  
        

        echo " "
        mkdir ${sampleName}
        #
        # find all the fastq files and make sure they are uncompressed
        #
        if [[ ${tarFile} == *".tar" || ${tarFile} == *".tar.gz" ]]; then
            tar -xvf ${tarFile} --directory=${sampleName}
            listOfFastQFiles=`find ${sampleName} -name "*fastq*"`

            
            printf "\nckecking if quant files need to be uncompressed. \nlistOfFastQFiles: $listOfFastQFiles\n"

            #
            # test if gz and uncompress
            #
            for quantFile in $listOfFastQFiles;
        do
        aedwip do not use gzip -t. it is slow see cut-n-paste.sh it uses file and grep
                gzip -t $quantFile 2>/dev/null
                if [ $? -eq 0 ];
                then
                    gzip -d $quantFile &
                    # else
                    #     not a compressed file
                fi

                printf "\n ****** debug\n"
                jobs
            done

            # wait for all background processes to complete
            # concurrent processing as possile
            wait
            
        fi

        printf "\n#########################\n"
        #
        # find the uncompress versions of all the fastq files
        #
        listOf_1_fastqFiles=`find ${sampleName} -name "*1.fastq*" `
        printf "\n listOf_1_fastqFiles\n $listOf_1_fastqFiles \n"

        listOf_2_fastqFiles=`find ${sampleName} -name "*2.fastq*" `        
        printf "\n listOf_2_fastqFiles\n $listOf_2_fastqFiles \n"


        #
        # concat the fastq files
        #
        cat $listOf_1_fastqFiles > ${sampleName}.1.fastq
        cat $listOf_2_fastqFiles > ${sampleName}.2.fastq

        #
        # compress
        #
        gzip ${sampleName}.1.fastq
        gzip ${sampleName}.2.fastq

    >>>

    output {
	# File firstEndFast = "${sampleName}.1.fastq"
	# File secondEndFast = "${sampleName}.2.fastq"
	File firstEndFast = "${sampleName}.1.fastq.gz"
	File secondEndFast = "${sampleName}.2.fastq.gz"
        
    }

    runtime {
	docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
	memory: "${memoryGb} GB"

        #
        # best practice 2 core min: one for os one for work
        # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#cpu
        # In Google Cloud: this is interpreted as "the minimum number of cores to use."
        #
	cpu: "2"
	disks: "local-disk ${diskSpaceGb} SSD"

        # if > 0 try and use preemptible cpus
        # integer value is number of times cromwell
        preemptible: '${runTimePreemptible}'        
    }

    meta {
        author: "Andrew E. Davidson"
	email: "aedavids@ucsc.edu"
    }    

}
