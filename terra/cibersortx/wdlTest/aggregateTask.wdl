version 1.0
workflow aggregateWorkflow {
    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

    call aggregateTask
}


task aggregateTask {
    #
    # CIBERSORTX mixture matrix assume each column is sample
    #
    # ref : https://xq-blog.dev/blog/wdl-for-loop/index.html
    #
    # generate a single data file. there will be n columns
    # The first column will be all 1. The second col will be all 2,
    # the nth col will be all n
    #
    # if we sum the columns will be easy to debug results
    #
    input {
        File csvDataFile
    
        # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/
        # String dockerImg = 'aedwip quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
        # docker pull alpine
        String dockerImg = 'aedavids/wdltest'

        String fileName = basename(csvDataFile)
    }

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
    echo "csvDataFile : ~{csvDataFile}"
    echo "fileName : ~{fileName}"
        
    # put copy of runtime parameters in output. Makes debug easier
    echo "runtime parameters"
    echo "dockerImg  : ~{dockerImg}"

    #
    # begin
    #
    set -euxo pipefail  # ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 

# weird indentation else we get a python indent error
sec=`python -c"""
import random
print( random.randint(10,30) )
"""`
        
    echo "AEDWIP sleep $sec"
        
    echo "!!!!!!!!!!! BEGIN aggregate task: sleep : $sec `date`"
    #outfile=`basename ~{csvDataFile}`
    python /bin/aggregate.py ~{csvDataFile} > '~{fileName}.aggregate.csv'

    echo "the PID of this process is $$"
    ls -l

    sleep $sec        
    echo "!!!!!!!!!!! END aggregate task: sleep: $sec `date`"             
    >>>

    output {
        #Array[File] aggregateCSV = glob("*.aggregate.csv")
        File aggregateCSV = '~{fileName}.aggregate.csv'
    }
    
    runtime {
        docker: '~{dockerImg}'
    }
    
    parameter_meta {
        #   https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#metadata-section
        csvDataFile: "A data file in comma separated format"
    }

    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

}

