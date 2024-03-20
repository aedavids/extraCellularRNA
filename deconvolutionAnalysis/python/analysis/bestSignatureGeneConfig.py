#
# bestSignatureGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

# BestSignatureGeneConfigCLI CommandLine display the doc string 
'''
The BestSignatureGeneConfig class should not be used from cli.
Use the cli for testing varg parsing. see main()
'''

from analysis.bestSignatureGeneConfigCLI import BestSignatureGeneConfigCommandLine 
import logging
import os
import pandas as pd
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2024-01-08'
__updated__ = '2024-01-08'

##############################################################################
class BestSignatureGeneConfig(SignatureGeneConfiguration):
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
        # use absolute value of log fold change to select best 
        # biologically signifigant genes
        # 
        significantDF = deseqDF.loc[selectSignificantRowsPS, :]
        absPS = significantDF['log2FoldChange'].abs()
        significantDF2 = significantDF.assign(absLog2FoldChange=absPS)    

        selectBestUpRegulatedRows = significantDF2.loc[:, 'absLog2FoldChange'] >= self.lfcThreshold
        significantDF3 = significantDF2.loc[selectBestUpRegulatedRows, :]
            
        significantDF3 = significantDF3.sort_values( by = ['baseMean'], ascending=False)

        self.logger.info(f'END')
        return significantDF3
    

################################################################################
def main(inCommandLineArgsList=None):
    '''
    use to test parsing of vargs
    '''
    # we only configure logging in main module
    loglevel = "WARN"
    loglevel = "INFO"
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    logger = logging.getLogger(os.path.basename(__file__))
    logger.warning("BEGIN")

    #
    # always log run time env to make debug easier
    #
    ORIG_PYTHONPATH = os.environ['PYTHONPATH']
    logger.warning(f'PYTHONPATH : {ORIG_PYTHONPATH}')
    logger.warning(f'FILE: {__file__}')
    logger.warning(f'PWD: {os.getcwd()}')

    cli = BestSignatureGeneConfigCommandLine(version=__version__ , 
                            author=__author__ ,
                            date=__date__, 
                            update=__updated__)
    
    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    design          = cli.args.design 
    padjThreshold   = cli.args.padjThreshold
    lfcThreshold    = cli.args.lfcThreshold
    dataSetName     = cli.args.dataSetName
    number          = cli.args.number
    title           = cli.args.title
    localCacheRoot  = cli.args.localCacheRoot


################################################################################
if __name__ == '__main__':
    '''
    hack used to test parsing of vargs
    '''
    main( inCommandLineArgsList=["--design", "tilda_gender_category",
                                 "--padjThreshold", "0.001",
                                 "--lfcThreshold", "2.0",
                                 "--dataSetName", "GTEx_TCGA",
                                 "--number" , "10",
                                 "--title", "a vs all",
                                 "--localCacheRoot", "myCacheRoot"
                                ])


