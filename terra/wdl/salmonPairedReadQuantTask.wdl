# ref:
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
# https://aedavids@github.com/aedavids/extraCellularRNA.git
# bin/runSalmon.pancreas.plasma.ev.long.RNA.sh
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
    String outDir

    String dockerImg = 'quay.io/biocontainers/salmon:0.14.1--h86b0361_1'
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

    
    command {
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

        #mkdir -p ${outDir} mkdir create a file that i can not remove when testing using cromwell

        # by convention foo.tar would have a root dir name foo. how ever we can not
        # guarantee conventions was followed
        # use sed remove the last slash
        # the not all docker images have GNU tar. works is use zcat
        refIndexDir=`zcat ${refIndexTarGz} | \
        tar -tf -  | \
        head -n 1 | \
        sed -e 's/\/$//'`
        
        # https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-mapping-based-mode
        # --libType A : automatically infer the library type
        # --gcBias: model will attempt to correct for biases in how likely a
        #      sequence is to be observed based on its internal GC content
        # --seqBias:  will attempt to correct for random hexamer priming bias, 
        #          which results in the preferential sequencing of fragments
        #          starting with certain nucleotide motifs.

        # AEDWIP  --recoverOrphans : only be used in conjunction with selective alignment)

        echo AEDWIP salmon quant \
        -i $refIndexDir \
        --libType A \
        -1 "${leftReads}" \
        -2 "${rightReads}" \
        -p 8 \
        --recoverOrphans \
        --validateMappings \
        --gcBias \
        --seqBias \
        --rangeFactorizationBins 4 \
        --output ${outDir}

        # should we gzip quant and tar aux_info? cmd_info.json can be helpful
    }

    output {
        # File cmd_info_json = '${sampleId}.cmd_info.json'
        # File quantFile     = '${sampleId}.quant.sf'

        #File aux_info     = '${sampleId}.aux_info.tar.gz'
        #File meta_json    = '${sampleId}.meta_info.json'
        #File resource_log = 'resource_usage.log'
        #File salmon_log   = '${sampleId}.salmon_quant.log'
    }

    runtime {
        # disks: 'local-disk ${diskSpaceGb} HDD'
        # cpu: '${runtime_cpu}'
        # memory: '${memoryGb} GB'
        docker: '${dockerImg}'

        # https://cloud.google.com/kubernetes-engine/docs/how-to/preemptible-vms
        # instances that last a maximum of 24 hours in general, and provide no availability guarantees.
        # Preemptible VMs are priced lower than standard Compute Engine
        # preemptible: '${runtime_preemptible}' 

    }
}
