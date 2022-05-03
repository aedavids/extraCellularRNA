#!/bin/sh

#
# return a tsv file. first col header name is 'Name' remain headers
# are the sample names
#

set -x 

for f in `ls *quant.sf*`;
do

    echo "\n****** $f"
    gzip -t $f 2>/dev/null
    if [ $? -eq 0 ];
    then
        echo uncompress
    else
        echo not a compressed file
    fi
    
done
