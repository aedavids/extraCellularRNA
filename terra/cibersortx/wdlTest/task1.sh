#!/bin/bash

n=10
set -x

# create a line of test data
# for i in $(seq 0 $n);
# do
#     echo $i >> line.txt
# done

# # create the test data file
# for i in $(seq 0 $n);
# do
#     cp line.txt $i;
# done

# listOfFiles=`ls |grep -v line.txt | sort -n`
# outFile="${n}.test.tsv"
# paste $listOfFiles > $outFile

for i in $(seq 0 $n);
do
    for j in $(seq 0 $n);
    do
        echo $i >> $i
    done
done

listOfFiles=`ls  | sort -n`
outFile="${n}.test.tsv"
paste $listOfFiles > $outFile
