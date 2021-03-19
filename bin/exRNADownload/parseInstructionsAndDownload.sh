#!/bin/bash

export scriptName=`basename $0`

if [ $# -ne 2 ]; then
    echo "error: usage $scriptName batchDownLoadInstructionFile rootDir"
    echo "example $scriptName testDownLoadInstr.aa ./kimLab/exRNA.org"
    echo "    testDownLoadInstr.aa is a batch of download parameters. see README.md for more info"
	echo "    all files and directory will be created under ./kimLab/exRNA.org"
	echo ""
    exit 1
fi

# track arguments and environment setting to make more reproduceable
echo cli: $0 $@

batchInstructionsFile="$1"
rootDir="$2"

input=$batchInstructionsFile 
while read -r line 
do 
	# strange parsing does not work correctly with tabs 
	line=`echo "${line}" | tr '\t' ','`
	
	# because the tokens are tab delimited and some the condition token 
	# may have spaces, we need to parse. If we do not parse the 
	# the downLoad script will get the wrong number of arguments 
 	downLoadFile=$(echo $line | cut -d "," -f 1) 
 	
 	url=$(echo $line | cut -d "," -f 2) 
 	
 	sampleId=$(echo $line | cut -d "," -f 3 ) 
 	
 	condition=$(echo $line | cut -d "," -f 4) 
 	# replace space with underscore
 	condition=`echo $condition | tr ' ' '_'` 
 	
 	assession=$(echo $line | cut -d "," -f 5 )  
 	
  	downLoadExRNA.orgData.sh ${rootDir} ${downLoadFile} ${url} ${sampleId} ${condition} ${assession} 
	exitStatus=$?
	if [ $exitStatus -ne 0  ]; then
		echo "ERROR: downLoadExRNA.orgData.sh ${rootDir} ${downLoadFile} ${url} ${sampleId} ${condition} ${assession}"
	fi
	
done < "$input" 
instrFile=`basename $batchInstructionsFile` 
dataIsUpSMS.sh $phoneNumber $scriptName ${instrFile} exit status $exitStatus 
