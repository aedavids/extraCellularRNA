version 1.0
workflow cibersortxFractionsWorkflow {
    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

    call cibersortxFractionsTask
}

task cibersortxFractionsTask {
    input {
        String username
        String token
        File mixture
        File sigmatrix
        Int perm = 100
        String label = "fractions"
        Boolean QN
        Boolean verbose
        String dockerImg = 'aedavids/cibersortx_fractions'

       String mixtureFileName = basename(mixture)        
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
    echo "username : ~{username}"
    echo "token    : ~{token}"
    echo "mixture  : ~{mixture}"
    echo "sigmatrix : ~{sigmatrix}"
    echo "perm     : ~{perm}"
    echo "label    : ~{label}"
    echo "QN       : ~{QN}:"
    echo "verbose  : ~{verbose}"
        
    # put copy of runtime parameters in output. Makes debug easier
    echo "runtime parameters"
    echo "dockerImg  : ~{dockerImg}"

    #
    # begin
    #
    set -euxo pipefail  # ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 

    # can we call the end default entry point?

        
    echo "AEDWIP pwd : `pwd`" | tee aedwip.out

    qnArg="FALSE"
    if [ ~{QN} == true ]; then
        qnArg="TRUE"
    fi

        
    verboseArg="FALSE"
    if [ ~{verbose} == true ]; then
        verboseArg="TRUE"
    fi

    # CIBERSORTx expects the docker to have mounted the
    # input and output directories to /src/outdir and /src/data
    mixtureFileName=`basename ~{mixture}`
    ln -s ~{mixture} /src/data/
        
    sigmatrixFileName=`basename ~{sigmatrix}`
    ln -s ~{sigmatrix} /src/data/

    ls -l /src/data
        
    /bin/sh -c "/src/CIBERSORTxFractions\
        --username ~{username} \
        --token ~{token} \
        --mixture $mixtureFileName \
        --sigmatrix $sigmatrixFileName \
        --perm ~{perm} \
        --label ~{label} \
        --QN $qnArg \
        --verbose $verboseArg"
        
   # cromwell expects output to be in current working directory
   ls /src/outdir
   cp -r /src/outdir/* .


   # we want to include the fileName in the output so that our merge task
   # can assemble the results in the correct order
   mv 'CIBERSORTx_~{label}_Results.txt' 'CIBERSORTx_~{label}_~{mixtureFileName}.Results.txt'
        
   ls .
    >>>

    output {
        # fractions output file will be of form 'CIBERSORTx_' + ~{label} + '_Results.txt'
        # example CIBERSORTx_aedwip_label_Results.txt
        #File fractionsResults='CIBERSORTx_~{label}_Results.txt'
        File fractionsResults='CIBERSORTx_~{label}_~{mixtureFileName}.Results.txt'


    }
    
    runtime {
        docker: '~{dockerImg}'
    }
    
    # parameter_meta {
    #     #   https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#metadata-section
    #     n: "size of the test data set"
    # }

    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

}

