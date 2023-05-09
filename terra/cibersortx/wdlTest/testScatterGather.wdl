version 1.0
#
# proof of concept.
# 1. generate a test file
# 2. split into parts
# 3. use scatter to process parts in parallel
# 4. user gather to combine the parts into a single results
#
# ref:
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
# https://wdl-docs.readthedocs.io/en/stable/WDL/ScatterGather_parallelism/
# https://xq-blog.dev/blog/wdl-for-loop/index.html


# use zip file to import local wdl files#
import "aggregateTask.wdl"  as aggTask
import "createTestDataTask.wdl"  as createTask
import "partitionDataTask.wdl" as partTask
import "mergeTask.wdl" as mergeTask

workflow testScatterGatherWorkflow {
    input {
        Int n
        Int numSamplesInPartition
    }

    call createTask.createTestFile {
        input:
        n=n
    }

    call partTask.partitionDataTask {
        input:
        #numSamplesInPartition=testScatterGatherWorkflow.numSamplesInPartition,
        numSamplesInPartition=numSamplesInPartition,
        #numSamplesInPartition=4,                
        dataFile=createTestFile.testData
    }

    #'scattter' ie call aggregate on parts in parallel
    #https://wdl-docs.readthedocs.io/en/stable/WDL/ScatterGather_parallelism/

    scatter (part in partitionDataTask.dataPartitions) {
        call aggTask.aggregateTask {
            input:
            csvDataFile=part
            }
    }

    # # aedwip agg outputs file wdl name aggregateCSV it is input to mergeTask.
    # # we do not explicitly use the gather keyword. It is implicit
    # call mergeTask.mergeTask {
    #     input:
    #     csvParts= partitionDataTask.dataPartitions
    # }

    # aedwip agg outputs file wdl name aggregateCSV it is input to mergeTask.
    # we do not explicitly use the gather keyword. It is implicit
    call mergeTask.mergeTask {
        input:
        parts=aggregateTask.aggregateCSV
    }


    output {
        # workflow output
        #String aedwipWorkFlowOut = "HELLO WORLD"
        File results = mergeTask.results
        }
    
    meta {
        author: "Andrew E. Davidson"
        email: "aedavids@ucsc.edu"
    }
    
    parameter_meta {
        n: "size of the test data set to generate"
        numSamplesInPartition:  "The maxium number of samples in each parition"
        }
}

    

