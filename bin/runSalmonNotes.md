# rough design for Salmon processing

1) assume we have have access to fastq files
cd $kl/pancreas.plasma.ev.long.RNA
rooDir=`pwd`

$ ll fastq | head
total 1.5T
-rw-r--r-- 1 aedavids kimlab  1.2G Sep 30 20:39 SRR10080507_pass_1.fastq.gz
-rw-r--r-- 1 aedavids kimlab  1.3G Sep 30 20:39 SRR10080507_pass_2.fastq.gz
-rw-r--r-- 1 aedavids kimlab  2.3G Oct  3 22:51 SRR10080508_pass_1.fastq.gz
-rw-r--r-- 1 aedavids kimlab  2.5G Oct  3 22:51 SRR10080508_pass_2.fastq.gz
-rw-r--r-- 1 aedavids kimlab  2.3G Oct  4 10:55 SRR10080509_pass_1.fastq.gz
-rw-r--r-- 1 aedavids kimlab  2.5G Oct  4 10:55 SRR10080509_pass_2.fastq.gz
-rw-r--r-- 1 aedavids kimlab  2.7G Oct  4 00:15 SRR10080510_pass_1.fastq.gz
-rw-r--r-- 1 aedavids kimlab  2.9G Oct  4 00:15 SRR10080510_pass_2.fastq.gz
-rw-r--r-- 1 aedavids kimlab  1.3G Oct  4 00:53 SRR10080511_pass_1.fastq.gz
(base) [aedavids@courtyard pancreas.plasma.ev.long.RNA]$ 


TODO: move each SRR*pass_*.fastq.gz into a separate dir. This allows us 
to group all the processing done at the single sample level. E.G. salmon

$ mkdir ${rootDir}/SRR10080507
$ mv ${rootDir/fastq/SRR10080507* ${rootDir}/SRR10080507/

do we want to use the docker image?
docker pull combinelab/salmon:latest

https://combine-lab.github.io/salmon/getting_started/

2) create the salmon index file for TE if it does not exist
   salmonIndex = AEDWIP # ????gen.32.ucsc.rmsk.index???
   salmonIndex.fa.gz = AEDWIP
   if ! /public/groups/kimlab/indexes/${salmonIndex}
	cd /public/groups/kimlab/indexes/
	salmon index -i  ${salmonIndex} -t ${salmonIndex.fa.gz}


3) for a each part file
       for each sample
       	   salmonOut=AEDWIP
       	   cd 
	   # https://salmon.readthedocs.io/en/latest/salmon.html#description-of-important-options
	   # -p 6 : 6 threads
	   salmon quant -i ${salmonIndex} -p 6 --libType A \
	   --gcBias --biasSpeedSamp 5 \
	   -1 sample_01_1.fastq.gz -2 sample_01_2.fastq.gz \
	   -o ${salmonOut}

!!!!!!!!!!!!! How do we get the meta data encoded in the 'name' column ?

ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|lncRNA|


how did this get incorporated into quant.sf

tx2Gene.csv <- file.path( '/home/kimlab/genomes.annotations/gencode.32', 
                          'gen.32.tx.2.gene.csv')
  gen.32.tx.2.geneDF <- read.table(file=tx2Gene.csv, header=FALSE, sep ="|")




#####################
cancer, health


mkdir SRR10080507/salmon.te/
mkdir SRR10080507.cancer/salmon.te/
mkdir SRR10080508.health/salmon.te/



/public/groups/kimlab/indexes/sel.aln.gen.34.ucsc.rmsk.index.salmon.1.2.1





2:29
/public/home/rreggiar/projects/aale.kras/scripts/{salmonRun.sh, trimmomaticRun.sh, starRun.sh}