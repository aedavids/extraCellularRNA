#!/bin/sh

#
# return a tsv file. first col header name is 'Name' remain headers
# are the sample names
#

if [ "$#" -lt 1 ]; then
    echo "ERROR arguments are not correct"
    echo "Usage: $0 outputFile"
    echo "create a tsv file from all the salmon quant.sf files in the current directory"
    echo "will uncompress if in gz format"
    echo " example $0 myMatrix; will produce a file myMatrix.tsv"
    exit 1
fi

outputFile=$1

# ref: https://gist.github.com/vncsna/64825d5609c146e80de8b1fd623011ca 
#set -euxo pipefail
set -x


# check if file is compressed or not
for f in `ls *quant.sf*`;
do

    printf "\n****** $f"
    gzip -t $f 2>/dev/null
    if [ $? -eq 0 ];
    then
        gzip -d $f
    else
        printf not a compressed file
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
mkdir -p cut.out
for q in `ls *.quant.sf`;
do
    printf "\n****** $q"
    sampleName=`echo $q | cut -d . -f 1`
    sed '1d' $q | cut -f 5 > cut.out/${sampleName}
done

#
# get the transcript names
#
quantFile=`ls *.quant.sf | head -n 1`
sed '1d' ${quantFile} | cut -f 1  > names.txt

#
# combine into a single table
#
paste names.txt `ls cut.out/* | sort` > tmpTable.tsv

#
# reconstruct the header line
#

printf "Name" > header.txt
for s in `ls cut.out/* | sort`;
do
    sampleName=`basename $s`
    printf "\t${sampleName}" >> header.txt
done
printf "\n" >> header.txt

# strip out all the new lines to convert column into a row
#sed -e 's/\n//g' column.txt > rowHeader.txt

cat header.txt tmpTable.tsv > "${outputFile}.tsv"
