#!/bin/bash

# aedavids@ucsc.edu
# ref: /public/home/rreggiar/projects/aale.kras/scripts/salmonRun.sh


function runSalmon() {
    # runs salmon on one sample and outputs to that directory
    inputDir="$1"
    salmonIndex="$2"
    outputDir="$3"
    outputPath="$inputDir"/"$outputDir"
    
    set -x # turn debug on
    # set +x # turn debug off

    if [[ ! -f "$outputPath"/quant.sf ]]; then

        mkdir "$outputPath"

        if [[ -f "$inputDir"/output_single_end.fq.gz ]]; then

            trim_read="$inputDir"/output_single_end.fq.gz

            salmon quant \
                -i "$salmonIndex" \
                --libType A \
                -r "$trim_read" \
                -p 8 \
                --validateMappings \
                --gcBias \
                --seqBias \
                --recoverOrphans \
                --rangeFactorizationBins 4 \
                --output "$outputPath" 

            exitStatus=$?
            if [ $? -ne 0 ]; then
                echo ERROR salmon "$inputDir" returned exit status "$exitStatus"
                continue
            fi

            #set -x # turn debug on
            set +x # turn debug off

        else

            trim_fwd="$inputDir"/output_forward_paired.fq.gz
            trim_rev="$inputDir"/output_reverse_paired.fq.gz


            salmon quant \
                -i "$salmonIndex" \
                --libType A \
                -1 "$trim_fwd" \
                -2 "$trim_rev" \
                -p 8 \
                --validateMappings \
                --gcBias \
                --seqBias \
                --recoverOrphans \
                --rangeFactorizationBins 4 \
                --output "$outputPath" 

            exitStatus=$?
            if [ $? -ne 0 ]; then

                echo ERROR salmon "$inputDir" returned exit status "$exitStatus"
                continue

            fi
        fi

    fi    
}

