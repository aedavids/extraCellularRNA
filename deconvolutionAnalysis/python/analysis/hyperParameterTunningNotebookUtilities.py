#
# hyperParameterTunningNotebookUtilities.py
# refactored some common jupyter notebook code
#
# Andrew E. Davidson
# aedavids@ucsc.edu
# 1/14/24
#

import ipynbname

# use display() to print an html version of a data frame
# useful if dataFrame output is not generated by last like of cell
from IPython.display import display

import numpy as np
import pandas as pd
# display all columns and rows
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

import pathlib as pl
import pprint as pp
import os
import sys

from analysis.findMisclassificationErrors import findMisclassificationErrors
from analysis.hyperParameterTunningMetrics import metricsRunner, elifeCols, lungCols
from analysis.hyperParameterTunningMetrics import findFile, findSummaryMetricsCols
from analysis.utilities import findAllCategories, findAllGenes
from analysis.utilities import findIntersectionsWithDegree
from analysis.utilities import loadDictionary, saveDictionary

root = "/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category"
notebookName = ipynbname.name()
outDir = f'{root}/hyperParameter/{notebookName}.out'
print( f'output dir: \n{outDir}' )
os.makedirs(outDir, exist_ok=True)

################################################################################
def evaluate(
    root: str,
    outDir: str,    
    resultsDirs : list[str],
    outFilePrefix : str,
    metric : str,
    stageName : str,
    threshold : float,
    verbose : bool = False,
    ) -> tuple[pd.DataFrame, pd.DataFrame] :

    retDF, retBellowThresholdDF = metricsRunner(root, outDir, outFilePrefix, resultsDirs, 
                        metric=metric, threshold=threshold, verbose=verbose)

    display( retDF.loc[:, findSummaryMetricsCols(metric) + elifeCols  ] )

    print(f'\n {stageName} classes < {threshold} {metric}')
    selectRowsBellow = retBellowThresholdDF.loc[:, "stage"] == stageName
    
    display( retBellowThresholdDF.loc[selectRowsBellow, 'category'] )

    return (retDF, retBellowThresholdDF)
