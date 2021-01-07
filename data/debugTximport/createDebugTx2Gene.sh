#!/bin/bash

scriptName=`basename $0`
if [ $# -lt 1 ]; then
    echo "error: usage $scriptName  list of trascript  ids"
    echo "see findTranscriptId.sh"
    exit 1
fi

# create tmp file that will automatically be deleted is script crash
# uses standard system programming trick
# open file discriptor then unlink the file
#
# https://unix.stackexchange.com/a/181938
tmpfile=$(mktemp /tmp/${scriptName}.XXXXXX)
exec 3>"$tmpfile"
rm "$tmpfile"

# echo $tmpfile


tx2gene="genomes.annotations/gencode.v35/gencode.v35.ucsc.rmsk.tx.to.gene.csv"

for i in $@
do
    grep $i $tx2gene >> $tmpfile
done

# output unique entries
sort $tmpfile | uniq

rm "$tmpfile"
