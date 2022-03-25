#!/bin/sh

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 downLoadFilesList.txt outputDir"
    echo "each line in the file is the gsc cloud url to download"
    echo "sample url"
    echo "gs://anvil-gtex-v8-hg38-edu-ucsc-kim-lab-spark/quant/sparkGTEx-train.out/parseSalmonR\
eadsSelectCountsCLI.out/onlyNumReads/GTEX-111CU-0826-SM-5EGIJ"
    exit 1
fi

# set -x # turn debug on
# set +x # turn debug off

downloadFileList=${1}
outputDir=${2}

mkdir -p "${outputDir}"

set -x # turn debug on
# set +x # turn debug off

# gsutil cp behaves like unix cp it does not perserver the path
for url in `cat ${downloadFileList}`;
do
    sampleName=`basename ${url}`

    # remove '/' from end of url if it exists
    cleanURL=`echo ${url} | sed "s/\/$//"`
    gsutil -m cp -r "${cleanURL}" ${outputDir}
    status=$?

    if [ $status -ne 0 ]; then
        echo "WARNING gsutil cp returned status $status for $sampleName"
    fi
    
        
    
    echo "downloaded $sampleName"    
done

echo "\n\nfinished"
