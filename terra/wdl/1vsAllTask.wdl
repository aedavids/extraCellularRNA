# ref:
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md

workflow deseq_one_vs_all {
    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }
    
    call one_vs_all
}


task one_vs_all {
    File countMatrix
    File colData
    String design
    String referenceLevel
    Boolean? isCSV

    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/
    # String dockerImg = 'aedwip quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
    String dockerImg = 'aedavids/edu_ucsc_kim_lab-1vsall_1.0'

    #
    # best practice 2 core min: one for OS one for work
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#cpu
    # In Google Cloud: this is interpreted as "the minimum number of cores to use."
    #
    # https://salmon.readthedocs.io/en/latest/salmon.html
    # We find that allocating 8 — 12 threads results in the maximum speed, threads
    # allocated above this limit will likely spend most of their time idle / sleeping
    #
    Int runTimeCpu = 8

    # random guess
    Int memoryGb = 64

    # random guess
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
    
    command <<<

    # we do not control the execution of this script
    # cromwell/docker ? puts set -x in stderr. this makes debugging harder
    # https://ops.tips/gists/redirect-all-outputs-of-a-bash-script-to-a-file/#using-exec-to-redirect-stdout-and-stderr
    # Redirect standard error to standard out such that 
    # standard error ends up going to wherever standard
    # out goes (the file).
    exec 2>&1

    # put copy on input parameters in output. Makes debug easier
    echo "input parameters"
    echo "countMatrix    : ${countMatrix} "
    echo "colData        : ${colData}"
    echo "design         : ${design}"
    echo "referenceLevel : ${referenceLevel}"
        
    # put copy of runtime parameters in output. Makes debug easier
    echo "runtime parameters"
    echo "memoryGb   : ${memoryGb}"
    echo "runTimeCpu: ${runTimeCpu}"
    echo "diskSpaceGb: ${diskSpaceGb}"
    echo "dockerImg  : ${dockerImg}"
    echo "runTimePreemptible: ${runTimePreemptible}"        


    set -x # turn shell trace debugging on 

        
    #
    # isCSV is a flag argument. --isCSV is a flag argument
    # flag or empty variable
    #
    isCSVTrue=${default="true" isCSV}
        isCSVFlag="--isCSV"

    # normally we want to use a file suffix to describe the
    # file format. Unfortunatly WDL does not allow us to
    # expand command variables in the output section
    # we need to know in advance what the output file name will
    # be
    # outFile="${referenceLevel}_vs_all.results.csv"
    outFile="${referenceLevel}_vs_all.results"
    if [ "$isCSVTrue" != "true" ]; then
        unset isCSVFLag
        # outFile="${referenceLevel}_vs_all.results.tsv"
    fi
        
    #
    # determin the concurrency level
    # do not create more threads 'BiocParallel child processes
    # than we have cores for. We need one core for OS
    #
    minRunTimeCPU=2
    if [  "${runTimeCpu}" -lt $minRunTimeCPU ]; then    
        echo "ERROR  ${runTimeCpu} must be >=  $minRunTimeCPU"
        exit 1
    fi
    numThr=$(expr "${runTimeCpu}" - 1)
        
    R CMD /home/rstudio/DESeqScript.R \
      --countMatrix "${countMatrix}" \
      --colData "${colData}" \
      --design "${design}" \
      --referenceLevel ${referenceLevel} \
      --outFile "$outFile" \
      --numCores $numThr \
      --estimateSizeFactorsOutfile "estimatedSizeFactors.csv" \
      --oneVsAll \
      $isCSVFlag 2>&1 

    >>>


    output {
        File outFile="${referenceLevel}_vs_all.results"
    }


    runtime {
        disks: 'local-disk ${diskSpaceGb} SSD'
        cpu: '${runTimeCpu}'
        memory: '${memoryGb} GB'
        docker: '${dockerImg}'

        #
        # enable Out of Memory Retry
        # https://support.terra.bio/hc/en-us/articles/4403215299355
        # https://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#maxretries
        # https://cromwell.readthedocs.io/en/develop/cromwell_features/RetryWithMoreMemory/#retry-with-more-memory
        #
        # maxRetries: 3         

        # if > 0 try and use preemptible cpus
        # integer value is number of times cromwell
        preemptible: '${runTimePreemptible}'
    }

    parameter_meta {
        #   https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#metadata-section
        groupbyCountMatrix: "first column name is 'geneId' remain col names are the sample names"
        colData: "a table of sample information"
        design: "indicates how to model the samples. example '~ age + sex + tissue_id' "
        referenceLevel: "tissue type present in the 'tissue_id' column in the colData table"
        isCSV: "boolean. default is true. if false files must be in TSV format"
    }

    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

}
