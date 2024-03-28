#
# bestLFCSignatureGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

# BestLFCSignatureGeneConfigCLI CommandLine display the doc string 
'''
The BestSignatureGeneConfig class should not be used from cli.
Use the cli for testing varg parsing. see main()
'''

#from analysis.bestSignatureGeneConfigCLI import BestSignatureGeneConfigCommandLine 
import logging
import os
import pandas as pd
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2024-03-22'
__updated__ = '2024-03-22'

##############################################################################
class BestLFCSignatureGeneConfig(SignatureGeneConfiguration):
    '''
    find best "n" genes. see findGenes() for details
    '''

    logger = logging.getLogger(__name__)

    ################################################################################
    def __init__(self, 
                    dataSetName : str, 
                    design : str, 
                    padjThreshold : float, 
                    lfcThreshold : float, 
                    n : int, 
                    localCacheRootPath : str, 
                    title : str)  :
        super().__init__(
                    dataSetName=dataSetName,
                    design=design,
                    padjThreshold=padjThreshold,
                    lfcThreshold=lfcThreshold,
                    n=n,
                    localCacheRootPath=localCacheRootPath,
                    title=title
        )
    
        '''
        see pipeline.dataFactory.signatureGeneConfig.SignatureGeneConfig().__init__ 
        for documentation
        '''

        self.logger.info(f'BEGIN')
        self.logger.info(f'END')

   ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        '''
        Find genes that that are statistically signifigant with  lfc <= -2.0 or >= 2.0
        and baseMean >= median( baseMean )
        sorted by lfc 

        This algorithm selects biomarkers with the strongest 'signal' not biomarkers that
        are the most differentially expressed

        arguments:
            deseqDF:
                results of DESeq2 as a pandas dataframe            

                
            fileName:
                useful for logging, debugging, and addition down stream processing   

        return:
            pandas dataframe
            top N biologically signifigant genes sorted by baseMean
            
        ref: findBestSignatureGenes(deseqDF, signatureGeneConfig)
            extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
            bestSignatureGeneConfig.py
            extraCellularRNA/deconvolutionAnalysis/jupyterNotebooks/explore1vsAllResults.ipynb
        '''
        self.logger.info("BEGIN")

        colsToReturn = deseqDF.columns
        sortedDF = self._select(deseqDF, fileName)

        # return the top n
        retDF = sortedDF.loc[:, colsToReturn].head(n=self.n)

        self.logger.info("END")
        return retDF


    ################################################################################
    def _select(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        '''
        select genes by padj and logfold changes and order by baseMean
        '''
        self.logger.info(f'BEGIN')
        
        # colsToReturn = deseqDF.columns
    
        #
        # find statistically signifigant genes
        #
        selectSignificantRowsPS = deseqDF.loc[:,"padj"] < self.padjThreshold
        self.logger.info("number of genes with padj < {} : {}".format(self.padjThreshold,
                                                        selectSignificantRowsPS.sum()))

        #
        # use absolute value of log fold change to select best up or down regulated
        # biologically signifigant genes
        # 
        significantDF = deseqDF.loc[selectSignificantRowsPS, :]
        absPS = significantDF['log2FoldChange'].abs()
        significantDF2 = significantDF.assign(absLog2FoldChange=absPS)    

        selectBestLFCRows = significantDF2.loc[:, 'absLog2FoldChange'] >= self.lfcThreshold
        significantDF3 = significantDF2.loc[selectBestLFCRows, :]

        # select candidate biomarkers that a strong signal but maybe not the strongest signal
        baseMeanMedian = significantDF3.loc[:, "baseMean"].median()
        self.logger.info(f'baseMeanMedian :{baseMeanMedian}')

        selectAboveMedianRows = significantDF3.loc[:, "baseMean"] >= baseMeanMedian
        retDF = significantDF3.loc[selectAboveMedianRows, :]
            
        # we want the most differential expressed biomarkers
        retDF = retDF.sort_values(by="log2FoldChange", ascending=False)

        self.logger.info(f'END')
        return retDF

