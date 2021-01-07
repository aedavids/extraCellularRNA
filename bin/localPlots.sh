#!/bin/sh
#
# Andrew E. Davidson, aedavids@ucsc.edu
# 6/16/2020
#

#day.5.de.seq.csv* day.7.de.seq.csv

for i in `ls data/kras.ipsc/exo.data/*.csv`;
#for i in `ls data/kras.ipsc/*.csv`;
do
	
	dataFile=`basename $i`
	outFile="img/tmp/${dataFile}.png"
	set -x # turn debug on	
	python src/bme263DataVis/volcanoQuartilePlot.py -i $i -o $outFile --title $dataFile
	set +x # turn debug off	
	echo
done

#python src/bme263DataVis/volcanoPlot.py \
#		-i ${inputFile} \
#		-o ${outputFile} 
#		