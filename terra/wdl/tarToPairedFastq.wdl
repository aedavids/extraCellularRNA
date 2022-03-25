
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
        exec 2>&1

        echo "input values"
        echo "tar: ${tarFile}"        
        echo "sampleName: ${sampleName}"

        set -x # turn shell trace debugging on

        # create a sub dir and extract tar
        echo " "
        mkdir ${sampleName}
        # tar -xvf ${tarFile} -C ${sampleName}
        tar -xvf ${tarFile} --directory=${sampleName}

        # make it easier to debug
        echo " "
        ls -l ${sampleName}

        echo " "
        pwd
        ls -l .

        error=0
        unknownNameFirstEndFastq=`find ${sampleName} -name "*1.fast*"`
        if [ -z "$unknownNameFirstEndFastq" ]
        then
            echo "ERROR firstEndfastq not found"
            error=1
        else
            # we need to create a known name
        # cromwell can not expect command variable in output section
        #
        # if suffix is '*.1.fast' salmon produces error
        # ERROR: file [/cromwell_root/fc-secure-689c9432-cc55-446b-b247-25666c8ac96f/66c81aeb-1444-44b9-ac6e-3a3fe1794873/quantify/b09734c2-5fb3-41fe-9537-893c9405d7f4/call-tarToFastqTask/ESCA-2H-A9GF-TP.2.fast] has extension .fast, which suggests it is neither a fasta nor a fastq file (or gzip compressed fasta/q).
# Is this file compressed in some other way?  If so, consider replacing: 

# /cromwell_root/fc-secure-689c9432-cc55-446b-b247-25666c8ac96f/66c81aeb-1444-44b9-ac6e-3a3fe1794873/quantify/b09734c2-5fb3-41fe-9537-893c9405d7f4/call-tarToFastqTask/ESCA-2H-A9GF-TP.2.fast

# with

# <(decompressor /cromwell_root/fc-secure-689c9432-cc55-446b-b247-25666c8ac96f/66c81aeb-1444-44b9-ac6e-3a3fe1794873/quantify/b09734c2-5fb3-41fe-9537-893c9405d7f4/call-tarToFastqTask/ESCA-2H-A9GF-TP.2.fast)
#
#        which will decompress the reads "on-the-fly"
        
            #firstEndFast="${sampleName}.1.fast"
            firstEndFast="${sampleName}.1.fastq"        
            mv $unknownNameFirstEndFastq $firstEndFast
            echo "firstEndFast: $firstEndFast"
        fi

        unknownNameSecondEndFastq=`find ${sampleName} -name "*2.fast*"`
        if [ -z "$unknownNameSecondEndFastq" ]
        then
            echo "ERROR secondfastq not found"
            error=1
        else
            # we need to create a known name
            # cromwell can not expect command variable in output section        
            # secondEndFast="${sampleName}.2.fast"
            secondEndFast="${sampleName}.2.fastq"
            mv $unknownNameSecondEndFastq $secondEndFast
            echo "secondEndFast: $secondEndFast"
        fi

        if [ "$error" -ne 0 ]
        then
            echo "error command returns -1"
            exit 1
        fi
    >>>

    output {
	File firstEndFast = "${sampleName}.1.fastq"
	File secondEndFast = "${sampleName}.2.fastq"
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
