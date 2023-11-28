#
# upsetPlotFactory.py
#   transforms DESeq results files into format upsetplot expects
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
#
# argument passing may be a little funny. I converted structure prog impl to obj
# tried to make minimal number of changes to reduce testing requirements
#

import logging
import pandas as pd
import pathlib as pl
import pipeline.dataFactory as df
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration
# from plots.DESeqSelect import DESeqSelect
# import plots.DESeqSelect
#import plots.DESeqSelect as dd


import upsetplot as upsp

###############################################################################
class UpsetPlotDataFactory( object ):
    '''
    transforms DESeq results files into format upsetplot expects 

    ref: extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
    
    public functions. See doc string for details
        __init__()
        createUpsetPlotDataFromSetDict()
    '''
    logger = logging.getLogger(__name__)
    
    ################################################################################    
    def __init__(self):
        pass

    ################################################################################    
    def createUpsetPlotDataFromSetDict(self, setDict :dict) -> tuple[pd.DataFrame, dict]:
        '''
        TODO

        arguments:
            setDict
                    key: candidate signature filename. 
                        Example 'Whole_Blood_vs_all.results'

                    value: pandas dataframe with DESeq results for selected rows
                    ie header columns name,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj

                    ref: pipeline.dataFactory.driver.runSelectGenesOfInterest

        returns tuple[pd.DataFrame, dict]
            pandas multilevel index dataframe in upsetplot expected format

            intersection dictionary
                for each class, dataset, tissueId
                    key = tissueId
                    value = list[ names] . e.g. the gene names
        '''

        # format data to make it easy to create a multi index dataframe
        geneSetsDict = dict()
        for key, df in setDict.items():
            tissueId = key.replace('_vs_all.results', '')
            geneSetsDict[tissueId] = df.loc[:,"name"].values

        retUpspDF = upsp.from_contents(geneSetsDict)
        
        return (retUpspDF, geneSetsDict)
