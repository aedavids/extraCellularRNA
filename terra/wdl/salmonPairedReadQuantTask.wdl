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

    String sampleId
    File refIndexTarGz
    File leftReads
    File rightReads
    String outDir = "salmon.out"

    String dockerImg = 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
    #String dockerImg =  'ubuntu:latest'
    Int runtime_cpu = 8
    Int memoryGb = 8
    Int diskSpaceGb = 40

    #parameter_meta {
        #    library: 'Salmon library type: https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype; by default, automatically infer'
    #}

    # scatter TODO I think terra will automatically spin up a docker for each file
    # so no need to run scatter
    call salmon_paired_reads {
        input:
        sampleId=sampleId,
        refIndexTarGz=refIndexTarGz,
        leftReads=leftReads,
        rightReads=rightReads,
        outDir=outDir,

        dockerImg=dockerImg,
        runtime_cpu=runtime_cpu,
        memoryGb=memoryGb,
        diskSpaceGb=diskSpaceGb
    }
}

task salmon_paired_reads {
    String sampleId
    File refIndexTarGz
    File leftReads
    File rightReads
    String outDir

    String dockerImg
    Int runtime_cpu
    Int memoryGb
    Int diskSpaceGb

    
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

        set -x

        # AEDWIP salmon runs out of memory. Terra seems to ignore runtime configuration
        # debug
        echo "AEDWIP salmon runs out of memory runtime parameter memoryGb: ${memoryGb}"
        # print docker memory stats
        cat /sys/fs/cgroup/memory/memory.stat

        # what version of link are we using
        cat /etc/os-release
        
        # extract the actual tar file
        time zcat ${refIndexTarGz} | tar -xf -
        # not all distributions support -z
        # tar -xzf ${refIndexTarGz}

        export refIndexDir=`ls -t | head -n 1 | sed -e 's/\/$//' `

        ls -l $refIndexDir

        # make sure the extracted tar file and anything else cromwell copied
        # into our local bucket will always be removed
        # unlink ${refIndexTarGz}
        # AEDWIP our image does not contain unlink
        #tmpFiles=`find $refIndexDir`
        #for i in $tmpFiles;
        #do
        #        unlink $i
        #done

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
        #time {

            # AEDWIP debug terra runtime parameters normally we would not use sh
#            sh -c '\
            salmon quant \
              -i $refIndexDir \
              --libType A \
              -1 "${rightReads}" \
              -2 "${leftReads}" \
              -p 8 \
              --recoverOrphans \
              --validateMappings \
              --gcBias \
              --seqBias \
              --rangeFactorizationBins 4 \
            --output ${outDir} 
            #' &

             salmonRet=$?
        #}

        # AEDWIP check runtime memory usage
        # for i in {1..50};
        # do
        #     echo "debug memory stats i:$i"
        #     cat /sys/fs/cgroup/memory/memory.stat
        #     sleep 60
        # done
        
        echo "AEDWIP in time salmonRet=$salmonRetXXXX";

        
        if [ $salmonRet -eq 0 ]; then
          gzip -c ${outDir}/quant.sf > ${sampleId}.quant.sf.gz;
          tar -c ${outDir}/aux_info/*.gz > ${sampleId}.aux_info.tar.gz;
        else
          echo "Salmon ERROR code $salmonRet";
        fi

        
        # clean up tmp files
        #rm -rf $refIndexDir
     >>>

     output {
         File quantFile     = '${sampleId}.quant.sf.gz'
         File aux_info     = '${sampleId}.aux_info.tar.gz'

         # we can not return a directory. We have two choices
         # 1) tar the directory and return the tar as type 'File'
         # 2) return an array of file names ie Array[File]
         # glob does not return files in sub directories
         # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#globs
         #   Array[File] output_bams = glob("*.bam")
         #Array[File] salmon_aedwip = glob("${outDir}/*")
     }

     runtime {
         disks: 'local-disk ${diskSpaceGb} HDD'
         cpu: '${runtime_cpu}'
         memory: '${memoryGb} GB'
         docker: '${dockerImg}'

         # https://cloud.google.com/kubernetes-engine/docs/how-to/preemptible-vms
         # instances that last a maximum of 24 hours in general, and provide no availability guarantees.
         # Preemptible VMs are priced lower than standard Compute Engine
         # preemptible: '${runtime_preemptible}' 

     }
 }
