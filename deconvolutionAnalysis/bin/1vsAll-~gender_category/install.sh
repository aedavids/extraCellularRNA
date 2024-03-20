#
# Andrew E. Davidson
# aedavids@ucsc.edu
# 12/28/23
#
# hacky little script to install a list of pipeline run scripts
#

set -x

#l="best100Enriched_6_Degree1_selectiveEnrich_LUAD_LUSC_12.sh best100Enriched_6_Degree1_selectiveEnrich_LUAD_LUSC_3.sh* best100Enriched_6_Degree1_selectiveEnrich_LUAD_LUSC_6.sh* best100Enriched_6_Degree1_selectiveEnrich_LUAD_LUSC_9.sh"

#l="rank5GTEx_TCGA.sh rank10GTEx_TCGA.sh rank15GTEx_TCGA.sh rank20GTEx_TCGA.sh"
# git add $l

#l="best100Removed_5_GTEx_TCGA.sh best20Removed_5_GTEx_TCGA.sh best25Removed_5_GTEx_TCGA.sh best30Removed_5_GTEx_TCGA.sh best50Removed_5_GTEx_TCGA.sh"

#l="best200Removed_5_GTEx_TCGA.sh best500Removed_5_GTEx_TCGA.sh"

#l="best200RemovedGTEx_TCGA.sh best500RemovedGTEx_TCGA.sh"

#l="best25EnrichedGTEx_TCGA.sh best30EnrichedGTEx_TCGA.sh  best50EnrichedGTEx_TCGA.sh  best100EnrichedGTEx_TCGA.sh "

#l="best20Enriched_6_GTEx_TCGA.sh best25Enriched_6_GTEx_TCGA.sh best30Enriched_6_GTEx_TCGA.sh best50Enriched_6_GTEx_TCGA.sh best100Enriched_6_GTEx_TCGA.sh"

# l="best100EnrichedGTEx_TCGA.sh best100Enriched_6_GTEx_TCGA.sh"

#l="best200Degree1GTEx_TCGA.sh best100Degree1GTEx_TCGA.sh best50Degree1GTEx_TCGA.sh best30Degree1GTEx_TCGA.sh best25Degree1GTEx_TCGA.sh best20Degree1GTEx_TCGA.sh"

#l="best500Degree1GTEx_TCGA.sh"

#l="best100Degree1GTEx_TCGA.sh best200Degree1GTEx_TCGA.sh best20Degree1GTEx_TCGA.sh best25Degree1GTEx_TCGA.sh best30Degree1GTEx_TCGA.sh best500Degree1GTEx_TCGA.sh best50Degree1GTEx_TCGA.sh"

#l="best20FindAllDegree1_wl5.sh"

# l="best500FindAllDegree1_wl500.sh"

#l="best1GTEx_TCGA.sh best2GTEx_TCGA.sh best3GTEx_TCGA.sh best5GTEx_TCGA.sh best10GTEx_TCGA.sh"

# l="best10GTEx_TCGA.sh"
#l='best1CuratedDegree1.sh best2CuratedDegree1.sh best3CuratedDegree1.sh best5CuratedDegree1.sh best10CuratedDegree1.sh'

# return findGenes() do not artifically enrich: empty df if no degree1 genes for category
# bug fix did not sort in descending order
#l='best1CuratedDegree1.sh best2CuratedDegree1.sh best3CuratedDegree1.sh best5CuratedDegree1.sh best10CuratedDegree1.sh'

#l="best3CuratedDegree1.sh"

#l="best10CuratedDegree1_ce467ff.degree1LUSC1.sh best10CuratedDegree1_ce467ff.degree1LUSC2.sh best10CuratedDegree1_ce467ff.degree1LUSC3.sh best10CuratedDegree1_ce467ff.degree1LUSC5.sh"
root=/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category

l="best10CuratedDegree1_ce467ff.degree1LUSC10.sh"

for i in $l;
do
    runName=`basename $i .sh`
    echo "runName : $runName"
    outDir="${root}/${runName}/training"
    mkdir -p $outDir
    cp $i $outDir
    cp ~/extraCellularRNA/deconvolutionAnalysis/bin/pipeline.sh $outDir
done
    
