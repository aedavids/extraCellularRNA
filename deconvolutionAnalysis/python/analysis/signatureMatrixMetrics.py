#
# signatureMatrixMetrics.py
# Andrew E. Davidson
# aedavids@ucsc.edu
# 1/26/24
#
'''
public functions:
'''
import pandas as pd

################################################################################
def checkForSparcity(
        XDF : pd.DataFrame,
        metaDataDF : pd.DataFrame
    ) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.Series, float) :
    '''

    arguments
        XDF :
            a count dataframe where row index are sample ids, and columns are genes

        metaDataDF:
            a deseq colData DataFrame. with a 'sample_id' column name

    returns tuple with elements
        noZeroDF, 
            boolean: true if value != 0
        
        rowSumDF
            for each sample id, returns the number values ! = 0
         
        rowSumSummaryDF
            mean, std, and median for each category
        
        geneCountSeries
            for each genes returns the number of values ! = 0
        
        sparcity:
            percent of values that are not zero
    '''
    noZeroDF = XDF.loc[:,:] != 0
    # print(f'XDF.shape {XDF.shape}')
    # print(f'noZeroDF.shape {noZeroDF.shape}')

    byRows = 1
    rowSumSeries = noZeroDF.sum(axis=byRows)
    rowSumSeries.name = "rowSum"
    # print(f'len(rowSumSeries) : {len(rowSumSeries)}')

    rowSumDF = pd.merge(rowSumSeries, metaDataDF, how='inner', left_index=True, right_on='sample_id')
    rowSumDF = rowSumDF.loc[:, ["sample_id", "category", "rowSum"]]
    rowSumSummaryDF = rowSumDF.groupby("category")["rowSum"].agg(['mean', 'std', 'median']) 

    geneCountSeries = noZeroDF.sum()

    sparcity = geneCountSeries.sum() / (noZeroDF.shape[0] * noZeroDF.shape[1]) 

    return (noZeroDF, rowSumDF, rowSumSummaryDF, geneCountSeries, sparcity)