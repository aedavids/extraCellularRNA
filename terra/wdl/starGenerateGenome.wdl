# ref:
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md

workflow starGenerateGenome {
    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }
    
    call starGenerateGenomeTask
}

task starGenerateGenomeTask {
    File referenceFasta
    Int sjdbOverhang
    File? sjdbGTF_File
    String outputIndexFileName

    Int memoryGb = 100
    Int diskSpaceGb = 1000
    Int numThreads = 10

    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible
    # Take an Int as a value that indicates the maximum number of times Cromwell should request a
    # preemptible machine for this task before defaulting back to a non-preemptible one.
    # default value: 0    
    Int runTimePreemptible = 3

    command <<<
        # put copy input parms values in output to make debug easier
        echo "referenceFasta    : ${referenceFasta}"
        echo "sjdbOverhang      : ${sjdbOverhang}"
        echo "outputIndexFileName: ${outputIndexFileName}"
        echo "sjdbOverhang      : ${sjdbOverhang}"
        
        # put copy of runtime parameters in output. Makes debug easier
        echo ""
        echo "runtime parameters"
        echo "memoryGb   : ${memoryGb}"
        echo "numTreads: ${numThreads}"
        echo "diskSpaceGb: ${diskSpaceGb}"
        echo "runTimePreemptible: ${runTimePreemptible}"        

        # we do not control the execution of this script
        # cromwell/docker ? puts set -x in stderr. this makes debugging harder
        # https://ops.tips/gists/redirect-all-outputs-of-a-bash-script-to-a-file/#using-exec-to-redirect-stdout-and-stderr
        # Redirect standard error to standard out such that 
        # standard error ends up going to wherever standard
        # out goes (the file).
        exec 2>&1

        set -euxo pipefail  # ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca

        if [ ! -z ${sjdbGTF_File} ]
        then
            file  ${sjdbGTF_File} | grep gzip 2>&1 > /dev/null
            if [ $? -eq 0 ];
            then
                gzip -d ${sjdbGTF_File}
                gtfFileName=`basename  ${sjdbGTF_File} .gtf`
            # else
                #     printf not a compressed file
                gtfFileName=${sjdbGTF_File}
            fi

            GTF_FILE="--sjdbGTFfile $gtfFileName "
        fi

        mkdir -p  ${outputIndexFileName}

        STAR --runMode genomeGenerate runThreadN ${numThreads} \
            --genomeDir  ${outputIndexFileName} \
            --genomeFastaFiles ${referenceFasta} \
            --sjdbOverhang ${sjdbOverhang} \
            $GTF_FILE

        tar -cvzf ${outputIndexFileName}.tar.gz  ${outputIndexFileName}

        
    >>>

    output {
        File starIdx="${outputIndexFileName}.tar.gz"
    }
    
    runtime {
        docker: "us.gcr.io/tag-public/neovax-tag-rnaseq:v1"
        memory: "${memoryGb}GB"
        disks: "local-disk ${diskSpaceGb} SSD"
        cpu: "${numThreads}"
        preemptible: "${runTimePreemptible}"
    }

    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }
    
}





