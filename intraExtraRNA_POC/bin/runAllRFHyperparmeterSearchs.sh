# run all data sets used in intraExtraRNA_POC/jupyterNotebooks/elife/elifeBinaryRandomForestResults.ipynb

set -x

outdir="/private/groups/kimlab/aedavids/elife/hyperparmeterTunning"
mkdir -p $outdir
randomForestHyperparmeterStomachSearch.sh $outdir
randomForestHyperparmeterLungSearch.sh $outdir
randomForestHyperparmeterLiverSearch.sh $outdir

randomForestHyperparmeterCOAD_Search.sh  $outdir
randomForestHyperparmeterColon_TransverseSearch.sh $outdir
randomForestHyperparmeterColon_SigmoidSearch.sh $outdir

randomForestHyperparmeterEsophagus_ESCASearch.sh $outdir 
randomForestHyperparmeterEsophagus_MucosaSearch.sh $outdir
randomForestHyperparmeterEsophagus_MuscularisSearch.sh $outdir

