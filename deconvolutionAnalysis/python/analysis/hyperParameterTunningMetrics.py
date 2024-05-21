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
    retList = metricCols + ["numGenes", "numTypes", "numDegree1", "numAboveThreshold"]

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

        # add metrics to valuesDF
        valuesDF.loc[:, [f'mean_{metrics}'] ] = rowAverage
        valuesDF.loc[:, [f'std_{metrics}'] ] = rowStd
        valuesDF.loc[:, [f'median_{metrics}'] ] = rowMedian
        valuesDF.loc[:, ['numGenes'] ] = numGenes
        valuesDF.loc[:, ['numTypes'] ] = numTypes   
        valuesDF.loc[:, ['numDegree1'] ] = len( degree1Sets ) 
        valuesDF.loc[:, ['numAboveThreshold'] ] = numAboveThreshold
        
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
