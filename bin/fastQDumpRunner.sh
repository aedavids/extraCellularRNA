#!/bin/bash

# aedavids@ucsc.edu
## download tools and install
# - [SRA toolkit home page](https://ncbi.github.io/sra-tools/)
# - [download tool for centos](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz)
# - [how to use sra tools](https://edwards.sdsu.edu/research/fastq-dump/)


scriptName=`basename $0`
if [ $# -lt 1 ]; then
    echo "error: usage $scriptName  list of SRA tool sample ids"
    echo "example $scriptName  SRR10080507 SRR10080508"
    exit 1
fi

# log CLI and tool version
echo $@
echo `fastq-dump --version`

set -x
sraDir=/public/groups/kimlab/NCBI/sra
listOfSRA_IDS=${@:1}

function runner() {
	listOfSampleId="$1"
	for i in $listOfSampleId
	do
	    set -x # turn trace debug on

	    # sratools slow, generate errors and are hard to use
	    # make reboust by allowing re-runs
	    if [ ! -f "fastq/${i}_pass_1.fastq.gz" ] ;
	    then

		# download sra file
		# do not call fastq-dump with out calling prefetch first. there is good
		# probablity that fastq-dump connection will time out
		prefetch $i
		exitStatus=$?
		if [ $? -ne 0 ];
		then
		    echo ERROR prefetch $i returned exit status $exitStatus
		    continue
		fi

		# convert to fastq.gz file format
		fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-e --clip ${sraDir}/${i}.sra
		exitStatus=$?
		if [ $? -ne 0 ];
		then
		    echo ERROR fastq-dump ${sraDir}/${i}.sra returned exit status $exitStatus
		fi

		set +x # turn trace debug off

	    fi
	done
}


runner "${listOfSRA_IDS}"

