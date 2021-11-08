# ref:
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
# https://portal.firecloud.org/?return=terra#methods/mxhe/salmon_quant_array/9

# using womtool to validate
# $ java -jar ../../java/bin/womtool-58.jar validate salmon.wdl

workflow salmon_quant {
    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

    call salmon_paired_reads
}

task salmon_paired_reads {
    String sampleId
    File refIndexTarGz
    File leftReads
    File rightReads
    String outDir = "salmon.out"

    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/
    String dockerImg = 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
    #String dockerImg =  'ubuntu:latest'

    #
    # best practice 2 core min: one for os one for work
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#cpu
    # In Google Cloud: this is interpreted as "the minimum number of cores to use."
    #
    # https://salmon.readthedocs.io/en/latest/salmon.html
    # We find that allocating 8 — 12 threads results in the maximum speed, threads
    # allocated above this limit will likely spend most of their time idle / sleeping
    #
    Int runTimeCpu = 8

    # 64 was min size needed to all GTEx cases
    Int memoryGb = 64

    # reference is about 20 GB
    # 80 was determined by trial and error
    Int diskSpaceGb = 80

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
    
    command  <<<
        
        # see bash man page "SHELL BUILTIN COMMANDS" for details
        # set -euxo pipefail
        # set -e Exit immediately if a pipeline see shell builtin command it is more complicated
        # aedwip_set -u Treat unset variables and parameters as  an  error
        # set -x turn debug trace on. output goes to stderr, normal output goes to stdout; merge these
        # set -o print value of current options

        
        # we do not control the execution of this script
        # cromwell/docker ? puts set -x in stderr. this makes debugging harder
        # https://ops.tips/gists/redirect-all-outputs-of-a-bash-script-to-a-file/#using-exec-to-redirect-stdout-and-stderr
        # Redirect standard error to standard out such that 
        # standard error ends up going to wherever standard
        # out goes (the file).
        exec 2>&1

        # put copy of runtime parameters in output. Makes debug easier
        echo "runtime parameters"
        echo "memoryGb   : ${memoryGb}"
        echo "runTimeCpu: ${runTimeCpu}"
        echo "diskSpaceGb: ${diskSpaceGb}"
        echo "dockerImg  : ${dockerImg}"
        echo "runTimePreemptible: ${runTimePreemptible}"


        set -x # turn shell trace debugging on 

        # use cpuinfo to debug preemptiple/cpu/threading quotas and performance
        # severa have lots of cpu's this generates a lot of debug info
        # cat /proc/cpuinfo
        
        # how much disk has been used
        du -sh .
        
        # print docker memory stats
        #cat /sys/fs/cgroup/memory/memory.stat

        # what version of link are we using
        cat /etc/os-release
        
        # extract the actual tar file
        time zcat ${refIndexTarGz} | tar -xf -
        # not all distributions support -z
        # tar -xzf ${refIndexTarGz}

        export refIndexDir=`ls -t | head -n 1 | sed -e 's/\/$//' `

        ls -l $refIndexDir

        # https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-mapping-based-mode
        # --libType A : automatically infer the library type
        # --gcBias: model will attempt to correct for biases in how likely a
        #      sequence is to be observed based on its internal GC content
        # --seqBias:  will attempt to correct for random hexamer priming bias, 
        #          which results in the preferential sequencing of fragments
        #          starting with certain nucleotide motifs.

        # AEDWIP  --recoverOrphans : only be used in conjunction with selective alignment)
        mkdir -p ${outDir}
        
        # grouping command in bash using () cause them to run in a sub shell
        # using {} cause them to execute in the current shell
        # https://www.gnu.org/software/bash/manual/html_node/Command-Grouping.html
        # man time: can configured to display memory metric
        #time - {

            # AEDWIP debug terra runtime parameters run salmon in background so we can
            # collect system performance metrics
            #sh -c '\

            #
            # determin the number of threads.
            # do not create more threads than we have cores for
            # we need one core for OS
            #
            minRunTimeCPU=2
            if [  "${runTimeCpu}" -lt $minRunTimeCPU ]; then    
                echo "ERROR  ${runTimeCpu} must be >=  $minRunTimeCPU"
                exit 1
            fi
            numThr=$(expr "${runTimeCpu}" - 1)
            
            salmon quant \
              -i $refIndexDir \
              --libType A \
              -1 "${rightReads}" \
              -2 "${leftReads}" \
              -p $numThr \
              --recoverOrphans \
              --validateMappings \
              --gcBias \
              --seqBias \
              --rangeFactorizationBins 4 \
              --writeUnmappedNames \
            --output ${outDir} #\
            #' &

             salmonRet=$?
        #}

        # AEDWIP check runtime memory usage
        # seems to crash after about 45 min
        # for i in {1..50};
        # do
        #     echo "debug memory stats i:$i"
        #     date
        #     du -sh .
        #     cat /sys/fs/cgroup/memory/memory.stat
        #     sleep 300
        # done
        
        echo "AEDWIP in time salmonRet=$salmonRetXXXX";

        
        if [ $salmonRet -eq 0 ]; then
          gzip -c ${outDir}/quant.sf > ${sampleId}.quant.sf.gz;
          tar -c ${outDir}/aux_info/*.gz > ${sampleId}.aux_info.tar.gz;
        else
          echo "Salmon ERROR code $salmonRet";
        fi

        # how much disk has been used
        du -sh .        
        
        # clean up tmp files
        #rm -rf $refIndexDir
     >>>

     output {
         # we can not return a directory. We have two choices
         # 1) tar the directory and return the tar as type 'File'
         # 2) return an array of file names ie Array[File]
         # glob does not return files in sub directories
         # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#globs
         #   Array[File] output_bams = glob("*.bam")
         #Array[File] salmon_aedwip = glob("${outDir}/*")
         
         File quantFile     = '${sampleId}.quant.sf.gz'
         File aux_info     = '${sampleId}.aux_info.tar.gz'
     }

     runtime {
         disks: 'local-disk ${diskSpaceGb} SSD'
         cpu: '${runTimeCpu}'
         memory: '${memoryGb} GB'
         docker: '${dockerImg}'

         # https://cloud.google.com/kubernetes-engine/docs/how-to/preemptible-vms
         # instances that last a maximum of 24 hours in general, and provide no availability guarantees.
         # Preemptible VMs are priced lower than standard Compute Engine
         #
         # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
         # Take an Int as a value that indicates the maximum number of times Cromwell should request a
         # preemptible machine for this task before defaulting back to a non-preemptible one.
         # default value: 0
         # With a value of 1, Cromwell will request a preemptible VM
         preemptible: '${runTimePreemptible}'

        #
        # enable Out of Memory Retry
        # https://support.terra.bio/hc/en-us/articles/4403215299355
        # https://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#maxretries
        # https://cromwell.readthedocs.io/en/develop/cromwell_features/RetryWithMoreMemory/#retry-with-more-memory
        #
        maxRetries: 3         
     }
 }
