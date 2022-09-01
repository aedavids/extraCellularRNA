#!/bin/sh -x

clear

rm -i upsetPlot.zip;

createUpsetPlotZip.sh /Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA upsetPlot.zip;

echo
echo

PYTHONPATH="${PYTHONPATH}:./tmp/python"
python tmp/python/plots/geneSignatureUpsetPlot.py 
