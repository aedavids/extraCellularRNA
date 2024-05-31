#
# hyperParameterTunningMetrics.py
# refactored some common jupyter notebook code
#
# Andrew E. Davidson
# aedavids@ucsc.edu
# 12/28/23
#

# 
# public functions
# findClassesBellowThreshold
# findSummaryMetricsCols
# metricsRunner

# import fnmatch
import pandas as pd
# import pprint as pp
import numpy as np
# import os

from analysis.utilities import findAllCategories
from analysis.utilities import findAllGenes
from analysis.utilities import findFile
from analysis.utilities import findSetsWithDegree
from analysis.utilities import loadDictionary

import logging
logger = logging.getLogger( __file__ )

lungCols = ['Lung', 'LUAD', 'LUSC']

#
# elife categories
#
# (extraCellularRNA) aedavids@mustard $ pwd
# /private/groups/kimlab/alex/data/elife
# (extraCellularRNA) aedavids@mustard $ cut -d , -f 16 elife_scaled_metaData_2023-05-18.csv | sort | uniq -c
#      53 "Colorectal Cancer"
#      31 "Esophagus Cancer"
#      43 "Healthy donor"
#      26 "Liver Cancer"
#      35 "Lung Cancer"
#      36 "Stomach Cancer"

# elifePluseGTexCols = lungCols + [ "Colon_Sigmoid", "Colon_Transverse", "COAD", "READ", 
#              "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "ESCA", 
#              "Liver", "LIHC",
#               "Stomach", "STAD", ]

elifeCols = [
    # corresponding TCGA classes
    'LUAD', 'LUSC', # Lung Cancer
    "COAD", "READ", # "Colorectal Cancer"
    "ESCA", # "Esophagus Cancer"
    "LIHC", # Liver Cancer
    "STAD", # Stomach Cancer
    "Whole_Blood", # ?? Healthy donor
]

elifePluseGTexCols = elifeCols + ["Colon_Sigmoid", "Colon_Transverse", 
                                    "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis",
                                    "Liver", 
                                    "Stomach", ]
################################################################################
def adjacentRowSort(
        df : pd.DataFrame,
        col1 : str,
        col2 : str,
        numericCol : str,
        verbose : bool=False):

    '''
    give a data frame with 2 key columns and a numeric column return a data frame
    where rows with the same key are adjacent to each other and the sum of the adjacent 
    numeric columns are sum

    example input DataFrame col1=trueCat, col2=predCat, sumCol="errorCount

    trueCat                             predCat	                            errorCount	
    Skin_Sun_Exposed_Lower_leg	        Skin_Not_Sun_Exposed_Suprapubic	    116
    Colon_Transverse	                Colon_Sigmoid	                    90
    Breast_Mammary_Tissue	            Adipose_Subcutaneous	            84
    Skin_Not_Sun_Exposed_Suprapubic	    Skin_Sun_Exposed_Lower_leg	        82
    Esophagus_Muscularis	            Esophagus_Gastroesophageal_Junction	67
    Esophagus_Gastroesophageal_Junction	Esophagus_Muscularis	            66

    return a DataFrame of the form

    trueCat	                                predCat	                                errorCount_x errorCount_y
    Brain_Nucleus_accumbens_basal_ganglia	Brain_Caudate_basal_ganglia	            10		    55
    Brain_Caudate_basal_ganglia	            Brain_Nucleus_accumbens_basal_ganglia   45		    55
    Brain_Putamen_basal_ganglia	            Brain_Caudate_basal_ganglia	            27		    41
    Brain_Caudate_basal_ganglia	            Brain_Putamen_basal_ganglia	            14		    41
    Brain_Frontal_Cortex_BA9	            Brain_Anterior_cingulate_cortex_BA24	3		    38
    Brain_Anterior_cingulate_cortex_BA24	Brain_Frontal_Cortex_BA9	            35		    38

    arguments 
        df: a dataFrame with at least 2 columns.

        col1, col2: the key column names  

        numericCol : the column of numeric values to sum

        verbose :
            if True will leave the "sortKey" column in the return dataFrame
            this can make debugging easier
    '''
    logger.info("BEGIN")

    uniqClassNamesSet = set(df[col1]).union( set(df[col2]) )
    classNames = sorted( list(uniqClassNamesSet) )
    # start count = 1, prevent errors. example 0 + 3 == 1 + 2
    keyToGroupDict = {key: idx + 1 for idx, key in enumerate( classNames )}


    # add a key columns
    df.loc[:, 'sortKey'] = df.apply(adjacentRowSortKey, axis=1, args=(col1,col2, keyToGroupDict) )

    # Sort the dataframe by the sort key
    sortedDF = df.sort_values(by='sortKey')
    logger.info(f'sortedDF:\n{sortedDF}')

    # select rows that have the same sortKey and are adjacent
    # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.shift.html
    maskSeries = sortedDF['sortKey'] == sortedDF['sortKey'].shift(-1)
    # print( maskSeries[0:5] )
    # print(sum(maskSeries))

    # return both rows
    adjDF = sortedDF.loc[ maskSeries | maskSeries.shift(1), : ] 
    logger.error(f'AEDWIP aedwip return not')
    logger.info(f'adjDF:\n{adjDF}')

    # group adj rows and sum the numeric columns group
    sumDF = adjDF.groupby(by="sortKey")[numericCol].sum().reset_index()
    logger.info(f'sumDF:\n{sumDF}')


    # left join or outer join
    # returns all records from the left table, along with any matching records from the right table.
    # If there are no matching records in the right table, the result will include NULL values 
    # in the right columns for the records in the left table. Records from the right table that are 
    # not in the left table will not be included in the result
    resultDF = adjDF.merge(sumDF, on='sortKey', how='left', suffixes=(f'_{col1}', f'_{col2}'))
    logger.info(f'resultDF:\n{resultDF}')

    leftCol = f'{numericCol}_{col1}'
    sortedResultsDF = resultDF.sort_values(by=[leftCol, 'sortKey'], ascending=False)    

    if not verbose:
        sortedResultsDF.drop(columns='sortKey', inplace=True)

    logger.info("END")
    return sortedResultsDF

################################################################################
def adjacentRowSortKey(
        row : pd.Series, 
        col1 : str, 
        col2 : str, 
        keyToGroupDict : dict ):
    
    keyA = row[col1]
    keyB = row[col2]

    # print(f'keyA : {keyA} keyB : {keyB}')
    
    # Use the minimum of the indices of key_a and key_b as the primary sort key
    # This groups rows based on the keys
    defaultValue = len(keyToGroupDict)
    valueA = keyToGroupDict.get(keyA, defaultValue)
    valueB = keyToGroupDict.get(keyB, defaultValue)
    # groupKey = valueA + valueB error 4 + 1 == 3 + 2
    sortedList = sorted( (valueA, valueB) )

    # use hash value might not be str
    groupKey = hash( tuple(sortedList) )
    #print(f'a: {keyA}, b: {keyB} valueA : {valueA} valueB : {valueB} ret : {groupKey}\n')
    
    return groupKey

################################################################################
def findSummaryMetricsCols(metric : str) -> list[str] :
    '''

    arguments:
        metric:
            example "sensitivity"

    returns a list of columns
        example
            ['mean_specificity', 'std_specificity',	'median_specificity', 
            "numGenes", "numTypes", "numDegree1", "numAboveThreshold"]
    '''

    # mean_sensitivity	std_sensitivity
    metricCols = [f'{m}_{metric}' for m in ['mean', 'std', 'median'] ]
    retList = metricCols + ["numGenes", "numTypes", "numDegree1", "numAboveThreshold", "percentAboveThreshold"]

    return retList

################################################################################
def metricsRunner(
                rootDir : str,
                outDir : str,
                outFilePrefix : str,
                resultsDirs : list[str],
                metric : str = 'sensitivity',
                threshold : float = 0.7,
                verbose : bool = False,
                ) -> tuple[pd.DataFrame, pd.DataFrame] :
    '''
    TODO
    '''
    df, bellowThresholdDF = _loadMetrics(rootDir, resultsDirs, metrics=metric, threshold=threshold, verbose=verbose)
    outFile = f'{outDir}/{outFilePrefix}.{metric}.{threshold}.csv'
    print(f'\nsaving : {outFile}' )
    df.to_csv(outFile)

    outFile = f'{outDir}/{outFilePrefix}.{metric}.bellow.{threshold}.csv'
    print(f'\nsaving : {outFile}' )
    bellowThresholdDF.to_csv(outFile)

    return (df, bellowThresholdDF)

#
# private functions
#

################################################################################
def _loadInterectionDict(path: str, verbose : bool = False) -> dict :
    '''
    TOOD
    '''
    # intersectionDictPath = ! find $path -name "*.intersection.dict"
    intersectionDictPath = findFile(path, "*.intersection.dict")
    intersectionDictPath = intersectionDictPath[0]
    if verbose:
        print( f'\nload\n{intersectionDictPath}' )     
            
    retDict = loadDictionary(intersectionDictPath)

    return retDict

################################################################################
def _loadMetrics(rootDirPath: str, 
                resultsDirs : list,
                metrics : str = 'sensitivity',
                threshold : float = 0.7,
                verbose : bool = False,
                ) -> tuple[pd.DataFrame, pd.DataFrame] :
    '''
    todo

    arguments:
        metrics : 
            precision,recall,f1-score,support,specificity,sensitivity,tp,fn,fp,tn
            default = 'sensitivity'
        
        threshold : 
            report number of classes/types >= threshold
            default = 0.7

    returns (retDF, bellowThresholdDF)
    '''
    retDF = pd.DataFrame()
    retBellowDF = pd.DataFrame()
    for fn in resultsDirs:
        path = f'{rootDirPath}/{fn}'
        if verbose :
            print(f'path : {path}')

        # metricsPath = ! find $path -name "metricsRounded.csv"
        metricsPaths = findFile(path, "metricsRounded.csv")

        metricsPath = metricsPaths[0]
        if verbose:
            print( f'\nload {fn} :\n{metricsPath}' )
        df = pd.read_csv(metricsPath)

        cols = df.loc[ :,  'id' ]
        dataSeries = df.loc[ :, metrics ] 
        dataNP = dataSeries.values
        data = np.reshape(dataNP, (1, len(cols)) )
        tmpDF = pd.DataFrame(data, columns=cols, index=[fn] )
        
        # select all columns except 
        selectCols = ~tmpDF.columns.isin(['accuracy', 'macro avg', 'weighted avg'])
        valuesDF = tmpDF.loc[:, selectCols]

        # calculate additional descriptive stats
        byCols = 1
        rowAverage = valuesDF.mean(axis=byCols)
        rowStd     = valuesDF.std(axis=byCols)
        rowMedian  = valuesDF.median(axis=byCols)

        # how many genes where used?
        # signatureMatricPath = ! find $path -name "signatureGenes.tsv"
        #bug signatureMatricPaths = findFile(path, "*.intersection.dict")

        # signatureMatricPath = signatureMatricPaths[0]
        # if verbose:
        #     print( f'{signatureMatricPath}' ) 

        # df = pd.read_csv( signatureMatricPath, sep="\t" ) 
        # numGenes = df.shape[0]

        # how many sets have unique genes
        numTypes = valuesDF.head().shape[1]
        intersectionDict = _loadInterectionDict(path, verbose)
        degree1Sets = findSetsWithDegree(intersectionDict, 1)

        genes = findAllGenes(intersectionDict)
        numGenes = len(genes)

        # which type do not have degree 1?
        allCategoriesSet = findAllCategories(intersectionDict)
        if verbose:
            print(f'\n{fn} types without degree 1 intersections: \n {allCategoriesSet - degree1Sets}')

        # how many classes have metric >= threshold?
        rowsLogicalDF = valuesDF.loc[:, :] >= threshold
        numAboveThreshold = rowsLogicalDF.sum(axis=byCols)

        bellowThresholdDF = findClassesBellowThreshold(valuesDF, threshold)
        #print(f'\nclasses < threshold:\n{bellowThresholdDF}')

        percentAboveThreshold = numAboveThreshold / numTypes

        # add metrics to valuesDF
        valuesDF.loc[:, [f'mean_{metrics}'] ] = rowAverage
        valuesDF.loc[:, [f'std_{metrics}'] ] = rowStd
        valuesDF.loc[:, [f'median_{metrics}'] ] = rowMedian
        valuesDF.loc[:, ['numGenes'] ] = numGenes
        valuesDF.loc[:, ['numTypes'] ] = numTypes   
        valuesDF.loc[:, ['numDegree1'] ] = len( degree1Sets ) 
        valuesDF.loc[:, ['numAboveThreshold'] ] = numAboveThreshold
        valuesDF.loc[:, ['percentAboveThreshold'] ] = percentAboveThreshold
        
        retDF = pd.concat( [retDF, valuesDF] )
        retBellowDF = pd.concat( [retBellowDF, bellowThresholdDF])

    return (retDF, retBellowDF)

################################################################################
def findClassesBellowThreshold(
                                valuesDF : pd.DataFrame,
                                threshold : float
                                ) -> pd.DataFrame:
    '''
    returns a data frame of form. These are runs, categories that are bellow threshoold 

        row_index column_name  value
    0          0           A   True
    1          1           B   True
    2          2           C   True
    '''
    
    # valuesDF is a shape (1,83)
    logicalDF = valuesDF.loc[:, :] < threshold
   
    # Use stack() to get row indices and column names for True values
    resultsDF = logicalDF.stack().reset_index()
    resultsDF.index.name = "id"
    resultsDF = resultsDF.rename(columns={'level_0':"stage", 'id':'category', 0:"value"})

    # Filter to keep only rows where the value is True
    resultsDF = resultsDF[resultsDF['value']]

    return resultsDF
