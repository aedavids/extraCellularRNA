version 1.0
workflow createTestWorkflow {
    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

    call createTestFile
}


task createTestFile {
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
        Int n
    }
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/
    # String dockerImg = 'aedwip quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
    # docker pull alpine
    String dockerImg = 'aedavids/wdltest'

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
    echo "n : ~{n}"
        
    # put copy of runtime parameters in output. Makes debug easier
    echo "runtime parameters"
    echo "dockerImg  : ~{dockerImg}"

    #
    # begin
    #
    set -euxo pipefail  # ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 

    python /bin/createTestData.py ~{n} > '~{n}.test.csv'

    echo "AEDWIP pwd : `pwd`"
    ls -l
    >>>

    output {
        File testData='~{n}.test.csv'
    }
    
    runtime {
        docker: '~{dockerImg}'
    }
    
    parameter_meta {
        #   https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#metadata-section
        n: "size of the test data set"
    }

    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

}

