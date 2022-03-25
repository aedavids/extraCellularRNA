#!/bin/sh

if [ "$#" -lt 1 ]; then
    echo "ERROR"
    echo "Usage: $0 quantFileList.txt"
    echo "\t each line quantFileList.txt is the path to a file containing 1"
    echo "\t column extracted from a salmon quant.sf file"
    echo "\t all output is written to stdout"
    exit 1
fi

quantFileList=$1

set -x # turn debug on
# set +x # turn debug off


# single sample version

# # transpose numReads columns and save to file
# OUTFILE=t.csv
# head /scratch/aedavids/parseSalmonReadsSelectCountsCLI.out/GTEX-ZZPU-2726-SM-5NQ8O/part-00000-\
    # eb54d3ea-cfdf-4c86-a25e-a7a1282cb96f-c000.csv  | cut -d , -f 2 | tr  \\n , | sed "s/.$/\n/" > $OUTFILE

for quantFile in `cat $quantFileList`;
do
    # sample file 
    # sparkGTEx-train.out/parseSalmonReadsSelectCountsCLI.out/onlyNumReads/GTEX-1PIEJ-2026-SM-E6CPA/part-00000-02065580-92ea-456c-b5c6-37d6ec2bcc09-c000.csv
    #
    # head output
    # GTEX-1PIEJ-2026-SM-E6CPA
    # 0.0
    # 0.0
    # 468.014
    # 0.0

    # use cut to select first column (older version had 2 column file format)
    # use tr to transpose the column into a row. replace new line with ,
    # use sed to replace the last character, ie ',' with new line
    cut -d , -f 1 ${quantFile} | tr  \\n , | sed "s/.$/\\n/"
done



