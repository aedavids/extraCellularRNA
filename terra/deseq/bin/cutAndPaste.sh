#!/bin/sh

#
# return a tsv file. first col header name is 'Name' remain headers
# are the sample names
#

set -x 

# check if file is compressed or not
for f in `ls *quant.sf*`;
do

    echo "\n****** $f"
    gzip -t $f 2>/dev/null
    if [ $? -eq 0 ];
    then
        gzip -d $f
    else
        echo not a compressed file
    fi
    
done

#
# example of quant file
#
# $ head -n 3 READ-AF-2689-NT.quant.sf 
# Name	Length	EffectiveLength	TPM	NumReads
# ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|	1657	1490.000	0.000000.000
# ENST00000450305.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000002844.2|DDX11L1-201|DDX11L1|632|transcribed_unprocessed_pseudogene|	632	466.0000.000000	0.000

# cut the numReads column and remove the header line
for q in `ls *.quant.sf`;
do
    printf "\n****** $q"
    sed '1d' $q | cut -f 5 | head
done
