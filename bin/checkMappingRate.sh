#
# Andrew E. Davidson, aedavids@ucsc.edu
# 5/5/21
#

#
# check quality of our salmon ref
# we expect that salmon will have fewer reads than the original bam created by
# https://github.com/broadinstitute/gtex-pipeline

set -x # turn debug on
# set + x # turn debug off

scratch=/scratch/aedavids
dataDir=~aedavids/extraCellularRNA/data/terra/AnVIL_GTEx_V8_hg38_edu_ucsc_kim_lab

# choose sample that we ran through the salmonQuantWorkflow
sampleId=GTEX-111CU-0526-SM-5EGHK

aux_info=gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-salmon_paired_reads/cacheCopy/GTEX-111CU-0526-SM-5EGHK.aux_info.tar.gz

bam_file=gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/GTEX-111CU-0526-SM-5EGHK.Aligned.sortedByCoord.out.patched.md.bam

bam_index=gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/GTEX-111CU-0526-SM-5EGHK.Aligned.sortedByCoord.out.patched.md.bam.bai

firstEndFastq=gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-bamToFastq/cacheCopy/GTEX-111CU-0526-SM-5EGHK.1.fastq.gz

participant=GTEX-111CU

quant=gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-salmon_paired_reads/cacheCopy/GTEX-111CU-0526-SM-5EGHK.quant.sf.gz

secondEndFastq=gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-bamToFastq/cacheCopy/GTEX-111CU-0526-SM-5EGHK.2.fastq.gz

tissueId=Pancreas

tissueSiteDetail=Pancreas

unpairedFastq=gs://fc-secure-f5aa8a37-78e5-45f6-9c59-c643016f7d97/3f3873ae-ab43-4aec-a8dc-f5eae2b8e44b/quantify/e521dbca-0912-4269-9d22-4f12a78dbd25/call-bamToFastq/cacheCopy/GTEX-111CU-0526-SM-5EGHK.unpaired.fastq.gz


#
# find unmapped reads in original STAR aligned bam
#
sampleDir=${dataDir}/${sampleId}
mkdir -p $sampleDir

# I think we can not access this outside of terra
# printf "\n\n ######### fetch bam\n"
# localBamFile=${sampleDir}/${sampleId}.bam
# if [ ! -f $localBamFile ]; then
#     gsutil copy $bam_file $localBamFile
# fi



#
# get salmon files
#
localQuant=${sampleDir}/${sampleId}.quant.sf
if [ ! -f $localQuant ]; then
    gsutil copy $quant $localQuant
fi

#
# get aux_info
#
printf "\n\n ######### aux_info\n"
localAux_infoGZ=${sampleDir}/aux_info.tar.gz
if [ ! -f $localAux_infoGZ ]; then
    gsutil copy $aux_info $localAux_infoGZ
fi

localAux_infoDir=${sampleDir}/aux_info
tar -xvf ${localAux_infoGZ} -C  ${localAux_infoDir}

#
# get the fastq files
#

printf "\n\n########### fastq\n"
localFirstEndFastqGZ=${sampleDir}/${sampleId}.1.fastq.gz
if [ ! -f $localFirstEndFastqGZ ]; then
    gsutil copy $firstEndFastq $localFirstEndFastqGZ
fi

localSecondEndFastqGZ=${sampleDir}/${sampleId}.2.fastq.gz
if [ ! -f $localSecondEndFastqGZ ]; then
    gsutil copy $secondEndFastq $localSecondEndFastqGZ
fi


printf "\n\n########### salmon\n"
source /private/home/aedavids/unmappedReadsAnalysis/bin/salmonUnmapped.sh
salmonIndexDir="${scratch}/sel.align.gencode.v35.ucsc.rmsk.salmon.v1.3.0.sidx"
salmonOut="${sampleDir}/salmon.unmmaped.out"

if [ ! -d $salmonOut ]; then
    runSalmon $salmonIndexDir $localFirstEndFastqGZ $localSecondEndFastqGZ $salmonOut
fi


printf "\n\n ########### salmon mapping rate \n"
mappingRate=`grep 'Mapping rate' ${salmonOut}/logs/salmon_quant.log | tr "%" " "`
printf "mappingRate = ${mappingRate} \n"
