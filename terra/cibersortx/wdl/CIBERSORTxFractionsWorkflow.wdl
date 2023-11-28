version 1.0

# use zip to import local wdl files https://cromwell.readthedocs.io/en/stable/Imports/
import "cibersortxFractionsTask.wdl" as cft
#import "../wdlTest/partitionDataTask.wdl" as partTask
import "partitionDataTask.wdl" as partTask
#import "../wdlTest/mergeTask.wdl" as mergeTask
import "mergeTask.wdl" as mergeTask

workflow CIBERSORTxFractionsWorkflow {
    input {
        Int numSamplesInPartition
        String username
        String token
        File mixture
        File sigmatrix
        Int perm = 100
        String label = "fractions"
        Boolean QN
        Boolean verbose
        Boolean isCSV
    }


    call partTask.partitionDataTask {
        input:
        numSamplesInPartition=numSamplesInPartition,
        dataFile=mixture,
        isCSV=isCSV
    }

    #'scattter' ie call aggregate on parts in parallel
    #https://wdl-docs.readthedocs.io/en/stable/WDL/ScatterGather_parallelism/
    scatter( mixture in partitionDataTask.dataPartitions) {
        call cft.cibersortxFractionsTask {
            input :
            username = username,
            token = token,
            mixture = mixture,
            sigmatrix = sigmatrix,
            perm = perm,
            label = label,
            QN = QN,
            verbose = verbose
        }
    }

    # aedwip agg outputs file wdl name aggregateCSV it is input to mergeTask.
    # we do not explicitly use the gather keyword. It is implicit
    call mergeTask.mergeTask {
        input:
        parts= cibersortxFractionsTask.fractionsResults
    }

    output {
        # workflow output
        #String aedwipWorkFlowOut = "HELLO WORLD"
        File fractions = mergeTask.results
        }    
}
