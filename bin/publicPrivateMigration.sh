#!/bin/bash

# move /public/groups/kimlab to /private/groups/kimlab
# aedavids@ucsc.edu
#
# ref:
# notes on managing bulk down load
#   extraCellularRNA/bin/exRNADownload/README.md
# migration notes
#   - andy's google driver mydrive/kimLab/PRISM_ADMIN 
#   https://drive.google.com/drive/folders/1_gSbnvVAa6y3PFHGxga63m65fztCyqk7?usp=sharing
#
#   PRISM_adminNotes.md 
#   https://drive.google.com/file/d/1DXRL6Ayu33aDD5ivrTPNNTVbu1PlOc6M/view?usp=sharing
#

set -x # turn debug trace on
#set +x # turn debug trace off

# use export else variable will not be defined by inline script embedded in nohup
export scriptName=`basename $0`
andySMS=6508622639@txt.att.net
export phoneNumber=$andySMS

#
# configure log file
#
export dateStamp=`dateStamp.sh`
export logDir="./logs"
mkdir -p ${logDir}
export scriptLog="${logDir}/$scriptName.${dateStamp}.out"


#
# configure rsync parameters
#

# small test
# export folder=allen.institute
# export folder=aedavids

#export folder=aale.kras
#export folder=aale.luad.exo
#done export=aedavids
#export folder=allen.institute
#export folder=bin
#export folder=cleanUpDeepVarient.sh.ou
#export folder=covid.2019
#export folder=DavidL
#export folder=exoRNA-biomarkers-panc
#export folder=exosome-KRAS-inhib
#export folder=exoTIC-biomarkers
#export folder=genomes.annotations
#export folder=hiv-smRNA-exoTIC
#export folder=incoming.plasma.data
#export folder=indexes
#export folder=kras.ipsc
#export folder=NCBI
#export folder=panc.plasma.2020
#export folder=pancreas.plasma.ev.long.RNA
#export folder=plasma.ex.RNA.Patel
#export folder=plasma.microbiome
#export folder=RNAEditingIndexer
#export folder=rstudio-server
#export folder=SalmonTools
#export folder=seqData
#export folder=toil.tcga
#export folder=vikas
#export folder=welcome-to-the-kim-lab.wiki


export srcRoot=/public/groups/kimlab
export distRoot=/private/groups/kimlab

export srcDir=${srcRoot}/${folder}
#export distDir=${distRoot}/${folder}
export distDir=${distRoot}

#                                                                                           
# create a session                                                                          
# this will ensure jobs continue to run even after we                                       
# log out. It will also make easier to kill all the children                                
# if we need to restart for some reason                                                     
# 

# run everything in one shot
setsid sh -c 'set -x; \
           rsync -avhz aedavids@plaza:${srcDir} ${distDir}; \
           dataIsUpSMS.sh $phoneNumber $scriptName ${folder} exitStatus $? ' > $scriptLog 2>&1 &


