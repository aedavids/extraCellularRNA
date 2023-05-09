# https://github.com/openwdl/wdl/blob/main/versions/draft-2/SPEC.md#scatter--gather
# https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#scatter--gather
workflow loopWorkflow {
    call loopTask
}

# wdl this is draft 2. draft 2 is prior to wdl v 1
task loopTask {
        Array[String] csvParts = ["a", "b", "c"]
        Int n = 123
    
    command <<<
        set -euxo pipefail

        exec 2>&1
        echo "shel: $SHELL"
        $SHELL -h

        echo "AEDWIP The value of N should be 123 n = ${n}"

        for part in ${sep=" " csvParts}
        do
            echo "AEDWIP $part AEDWIP"
        done

        echo "\n\n********* create a bashArray"
        bashArray=('~{sep="' '" csvParts}')
        for part in "${bashArray[@]}"
        do
            echo "bashArray loop aedwip $part aedwip"
        done
        

        ls -l >results.csv
    >>>

    output {
        File  results_csv='results.csv'
    }
}


