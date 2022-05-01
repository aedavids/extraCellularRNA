#!/bin/bash


# get setNames
ls  signatureGenesValidate*.csv | grep -v lfcShrink | cut -d - -f 2 > setNames.txt
ls  signatureGenesValidate*.lfcShrink.csv | cut -d - -f 2 >  t
for i in `cat t`; do echo ${i}.lfcShrink >> setNamesLfcShrink.txt; done

# get num header lines
'rm' numHeaderLines.txt
for i in `cat setNames.txt`;
do
    echo 8 >> numHeaderLines.txt
done

'rm' numHeaderLines.LfcShrink.txt
for i in `cat setNamesLfcShrink.txt`;
do
    echo 7 >> numHeaderLines.LfcShrink.txt
done

         
# get file paths
'rm' files.txt
for i in `ls signatureGenesValidate*.csv | grep -v lfcShrink`;
do
    echo `pwd`/$i >> files.txt;
done


'rm' fileslfcShrink.txt
for i in `ls signatureGenesValidate*.lfcShrink.csv`;
do
    echo `pwd`/$i >> filesLfcShrink.txt;
done

# create csv file
paste -d , setNames.txt numHeaderLines.txt files.txt  > tmpDataSets.csv
paste -d , setNamesLfcShrink.txt numHeaderLines.LfcShrink.txt fileslfcShrink.txt  > tmpDataSetsLfcShrink.csv

echo "setName,numHeaderLines,filePath" > dataSets.csv
cat tmpDataSets.csv tmpDataSetsLfcShrink.csv >> dataSets.csv

'rm'  setNames.txt setNamesLfcShrink.txt files.txt filesLfcShrink.txt  tmpDataSets.csv tmpDataSetsLfcShrink.csv numHeaderLines.txt numHeaderLines.LfcShrink.txt
