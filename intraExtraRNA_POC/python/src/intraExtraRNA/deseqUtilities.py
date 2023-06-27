'''
Created on May 31, 2023

@author: andrewdavidson
aedavids@ucsc.edu

public functions:

def loadDESEqResults(
            filePath: pl.Path,
            skiprows: int = 0,
            index_col: str ='name'
        ) ->  pd.DataFrame :
        
def normalize(countDF:         pd.DataFrame, 
              scalingFactorDF: pd.DataFrame)-> pd.DataFrame:
                      
def selectBest(deseqDF:       pd.DataFrame, 
               topN:          int = 0, 
               padjThreshold: float = 0.001, 
               lfcThreshold:  float = 2.0) -> pd.DataFrame:
               
def selectGenesOfInterest(
            resultsDF : pd.DataFrame,
            geneList , # : list[str] requires python 3.9
            ) -> pd.DataFrame :

'''

import pathlib as pl
import pandas as pd

################################################################################
def loadDESEqResults(
            filePath: pl.Path,
            skiprows: int = 0,
            index_col: str ='name'
        ) ->  pd.DataFrame :
    '''
    reads filePath and filters out rows with baseMean <= 0
    
    arguments:
        filePath: path to DESEq output file
        
    skiprows
        int: default 0
        DESEq has multiple header lines. number of header lines to skip 
      
    index_col:
        string: default 'name'
        use this column in the DESEq results file as the data frame index
        
    return    
        pandas dataframe
    
    '''
    resultsDF = pd.read_csv(filePath,skiprows=skiprows, index_col=index_col)
    selectNonZeroBaseMeanRows =  resultsDF['baseMean'] > 0
    resultsDF = resultsDF.loc[selectNonZeroBaseMeanRows,:]
    
    return resultsDF

################################################################################
def normalize(countDF:         pd.DataFrame, 
              scalingFactorDF: pd.DataFrame)-> pd.DataFrame:
    '''
    countDF:
        type: pandas dataframe
        shape: number of genes x number of sample
        
    scalingFactorDF:
        type: pandas dataframe
        shape numSamples x 1
        values are floats
    '''
    # transpose shape number of samples x number of genes
    transposeGroupByDF = countDF.transpose(copy=True)
    
    # normalize counts
    # element wise multiplication . use values to to multiply a vector
    # print("before normalization")
    # display(transposeGroupByDF.head())
    
    normalizedDF = transposeGroupByDF * scalingFactorDF.values
    
    # print("after normalization")
    # display(normalizedDF.head())
    return normalizedDF

################################################################################
def selectBest(deseqDF:       pd.DataFrame, 
               topN:          int = 0, 
               padjThreshold: float = 0.001, 
               lfcThreshold:  float = 2.0) -> pd.DataFrame:
    '''
    1. select genes that that are statistically signifigant 
        adjusted p-value < padjThreshold   
        
    2. select genes that are biologically signifgant 
        lfc <= -lfcThreshold or >= lfcThreshold
        
    3. sort in descending order by baseMean. Hi baseMeans implies a stronger signal
       making them a better discriminator, and potentially more biologically 
       signifigant gene
    
    4. return the topN. If topN = 0 will return all
    
    arguments:
        deseqDF:
            results of DESeq2 as a pandas dataframe
            
        topN:
            integer. 
            default value = 0. 
            if topN is zero returns all
            
        padjThreshold:
            float
            default=0.01
            
        lcfThreshold:
            float
            default = 2.0
    
    return:
        pandas Dataframe
    '''    
    #
    # select statistically signifigant genes
    #
    selectSignificantRowsPS = deseqDF.loc[:,"padj"] < padjThreshold
    # print("number of genes with padj < {} : {}".format(padjThreshold,
    #                                                    selectSignificantRowsPS.sum()))
    significantDF = deseqDF.loc[selectSignificantRowsPS, :]    
    print(f'significantDF.shape : {significantDF.shape}')

    #
    # use absolute value of log fold change to select  
    # biologically signifigant genes
    # 
    absPS = significantDF['log2FoldChange'].abs()
    significantDF2 = significantDF.assign(absLog2FoldChange=absPS)    

    selectBestUpRegulatedRows = significantDF2.loc[:, 'absLog2FoldChange'] >= lfcThreshold
    significantDF3 = significantDF2.loc[selectBestUpRegulatedRows, :]
    # print("number of genes with absLog2FoldChange > {} : {}".format(lfcThreshold,
    #                                                    selectBestUpRegulatedRows.sum()))        
    #
    # genes with large baseMean have a strong signal, should be better discrimantors
    # and are more likely to be biologically signifigant
    #
    significantDF3 = significantDF3.sort_values( by = ['baseMean'], ascending=False)
    
    ret = significantDF3
    if topN > 0:
        ret = significantDF3.iloc[0:topN, :]

    return ret

################################################################################
def selectGenesOfInterest(
            resultsDF : pd.DataFrame,
            geneList , # : list[str] requires python 3.9
            ) -> pd.DataFrame :
    '''
    arguments:
        resultsDF : pandas data frame returned by loadDESEqResults()
        
        geneList: a list of genes name, E.G. T.E. names 
    '''
    resultNameList = resultsDF.index.to_list()
    resultNameSet = set(resultNameList)

    intersectionSet = resultNameSet.intersection(geneList)
    #print(len(intersectionSet))
    resultIdList = list(intersectionSet)
    
    retDF = resultsDF.loc[resultIdList, :]
    
    return retDF

################################################################################
if __name__ == '__main__':
    pass