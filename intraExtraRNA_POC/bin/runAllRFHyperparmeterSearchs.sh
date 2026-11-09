# run all data sets used in intraExtraRNA_POC/jupyterNotebooks/elife/elifeBinaryRandomForestResults.ipynb

set -x

outdir="/private/groups/kimlab/aedavids/elife/hyperparmeterTunning"
mkdir -p $outdir

# intraExtraRNA_POC/jupyterNotebooks/elife/elifeBinaryRandomForestResults.ipynb feature importance plots
# do not always have 10 genes
# search log files for
# logger.warning(f'Are all gene loci present XNP.shape : {XNP.shape}')
# if isinstance(df, pd.DataFrame):
#     logger.warning(f'Dataframe XNP.columns : {XNP.columns}')
 

# ran first 1/2 11/8/26 @ 12:31 pm
# randomForestHyperparmeterStomachSearch.sh $outdir
# randomForestHyperparmeterSTAD_Search.sh $outdir

# randomForestHyperparmeterLungSearch.sh $outdir
# randomForestHyperparmeterLUSCSearch.sh $outdir
# randomForestHyperparmeterLUADSearch.sh $outdir

# randomForestHyperparmeterLiverSearch.sh $outdir
# randomForestHyperparmeterLIHC_Search.sh $outdir

# ran second half 11/8/26 3:58
randomForestHyperparmeterCOAD_Search.sh  $outdir
randomForestHyperparmeterColon_TransverseSearch.sh $outdir
randomForestHyperparmeterColon_SigmoidSearch.sh $outdir

randomForestHyperparmeterEsophagus_ESCASearch.sh $outdir 
randomForestHyperparmeterEsophagus_MucosaSearch.sh $outdir
randomForestHyperparmeterEsophagus_MuscularisSearch.sh $outdir

