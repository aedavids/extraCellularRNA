version 1.0
workflow mergeWorkflow {
    input {
        Array[File] parts
    }
    
    call mergeTask {
        input:
        parts = parts
    }

    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }    
}


task mergeTask {
    input {
        Int n = 123
        Array[File] parts
    }
    
    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/
    # String dockerImg = 'aedwip quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
    # docker pull alpine
    String dockerImg = 'aedavids/wdltest'
    #String dockerImg = 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'

    
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
        
    # put copy of runtime parameters in output. Makes debug easier
    echo "runtime parameters"
    echo "dockerImg  : ~{dockerImg}"
        
    #
    # begin
    #
    set -euxo pipefail  # ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 

    #
    # wdl version 1.1 syntax
    # https://xq-blog.dev/blog/wdl-for-loop/index.html
    # arrays in wdl are not bash arrays
        # think of the bash portion of our task as a template script

    # we need to make sure the part are assembled in the correct order
    # previous tasks split the original sample into separate file
    # the file name can be sorted
    # we use unix file sort
    # step 1, create file that is a list of unsorted parts    
    for part in ~{sep=' ' parts}
    do
        echo "AEDWIP $part AEDWIP"
        echo "$part" >> parts.txt
    done

    # step 2 sort the parts and write to a new file in sort order
    python /bin/sortPartsFilePaths.py parts.txt > sortedParts.txt

   # combine the parts in the original order
   #cat `cat sortedParts.txt` > 'results.csv'

   # write the first part including its header
   cat `head -n 1 sortedParts.txt` > results.txt

   # strip the header line from the remaining files
   for partFile in `tail -n +2 sortedParts.txt`
   do
        tail -n +2 $partFile >> results.txt
   done
        
   >>>

    output {
        File  results='results.txt'
    }

    runtime {
        docker: '~{dockerImg}'
    }    

    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }

}
