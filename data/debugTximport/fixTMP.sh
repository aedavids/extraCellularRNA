#!/bin/bash

# turns out the TPM values are all zero
# quick hack to fix that 
#(base) [aedavids@plaza debugTximport]$ cut -f 4 ./pancreas.plasma.ev.long.RNA/data/healthy/SRR10080507/salmon.out/debugQuant.sf
#TPM
#0.000000
#0.000000
#0.000000
#0.000000



TMP=123.1234
files=`find . -name debugQuant.sf.save`
for i in $files
do 
    # we need to match on regex starting with tab
    # we need to fix the TPM and NumReads, TMP alone did not work,
    d=`dirname $i`
    #    sed s/\\t0.000000/\\t${TMP}/  $i #> $d/debugQuant.sf
    sed -e "s/\\t0.000000/\\t${TMP}/" -e  "s/\\t0.000/\\t55.0/" $i > $d/debugQuant.sf

    # I do not think we need to use BC
    #TMP=`echo $TMP + 100|bc`
done


