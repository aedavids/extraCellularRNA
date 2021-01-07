#!/bin/bash

# 
# find the debugQuant.sf files and select the name values
# this can be used to create debugVersion of tx2gene.csv
#

#set -x # turn debug on
# set + x # turn debug off


#
# the debug salmon quant.sf files only have 4 transcripts
# 2 start wtih ENST, the other are our custom annotations
# ex. hg38_rmsk_ALR/Alpha_range=chr22_KI270739v1_random:50071-73985_5'pad=0_3'pad=0_strand=-_repeatMasking=none
#

files=`find . -name "debug*.sf"`

for i in $files
do
    # use tail to select all rows except the header
    # The default delimiter for cut is tab
    tail -n 4 $i | cut -f 1  
done


