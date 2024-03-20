set -x
x="/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category"

#runNamesAll=`ssh mustard ls -d  "/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/{best,rank}*"`

runNamesRank=`ssh mustard ls -d  "/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/rank*"`
#for rn in `cat runNamesAll`;
for rn in $runNamesRank;          
do
    ldir=`basename $rn`
    mkdir $ldir;
    p="${rn}/training/*.out/upsetPlot.out"
    echo "p: ${p}"
    #ssh mustard "ls -d ${p}";
    scp -r  "mustard:${p}" $ldir
    echo ""

    p="${rn}/training/*.out/metrics"
    echo "p: ${p}"
    #ssh mustard "ls -d ${p}";
    scp -r  "mustard:${p}" $ldir


    p="${rn}/training/*.out/lungExploreClassificationErrors.sh.out"
    echo "p: ${p}"
    #ssh mustard "ls -d ${p}";
    scp -r  "mustard:${p}" $ldir

done

