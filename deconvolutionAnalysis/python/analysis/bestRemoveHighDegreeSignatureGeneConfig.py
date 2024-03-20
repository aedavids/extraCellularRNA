#
# bestRemoveHighDegreeSignatureGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

# BestRemoveHighDegreeSignatureGeneConfigCLI CommandLine display the doc string 
'''
The BestRemoveHighDegreeSignatureGeneConfig class should not be used from cli.
Use the cli for testing varg parsing. see main()
'''

import ast
import logging
import os
import pandas as pd

from analysis.bestRemoveHighDegreeSignatureGeneConfigCLI import BestRemoveHighDegreeSignatureGeneConfigCLI
from analysis.bestSignatureGeneConfig import BestSignatureGeneConfig
from analysis.utilities import loadDictionary
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-12-07'
__updated__ = '2023-12-07'

##############################################################################
class BestRemoveHighDegreeSignatureGeneConfig(BestSignatureGeneConfig):
    '''
    genes that are shared between many categories do not make good discrimnators

    find best "n" genes. see BestSignatureGeneConfig.findGenes() for details
    then remove genes that are elements of intersections with high degree.
    The degree of a set is the number of set that have elements in the intersection

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
                    title : str,
                    intersectionDictionaryPath : str,
                    degreeThreshold : int)  :
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

        additional arguments:
            intersectionDictionaryPath : str
                output from upstream stage of pipeilne. Dictionary key is multi-index identifying sets that share elements.
                value is the list of shared elements

                example: ./analysis/test/data/intersection.dict

                see plots.test.testUpsetPlots testIntersection()

            degreeThreshold : int
                the degree of an intersection is the number of sets that share elements
                genes from sets with degree > degreeThreshold will removed
        '''
        self.logger.info(f'BEGIN')

        self.intersectionDictionaryPath = intersectionDictionaryPath
        self.degreeThreshold = degreeThreshold

        # with open(self.intersectionDictionaryPath) as f: 
        #     data = f.read() 

        # # we can not convert the intersectionData to a DataFrame
        # # the size of the interesections is not constant. 
        # # we would need to pad the dictionary values
        # self.intersectionDict = ast.literal_eval(data)

        self.intersectionDict = loadDictionary(self.intersectionDictionaryPath)
        self.poorDiscriminators = self._findPoorDiscriminators()

        self.logger.info(f'END')

    ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        '''
        Find genes that that are statistically signifigant with  lfc <= threshold or >= threshold
        then removes genes from sets with high degree. These genes shared between many 
        categories and are not good descriminators.

        check log files to see which genes have been removed. They are probably 
        biologically interseting.
    
        arguments:
            deseqDF:
                results of DESeq2 as a pandas dataframe  

            fileName:
                useful for logging, debugging, and addition down stream processing   
    
        return:
            pandas dataframe
            
        ref: 
            BestSignatureGeneConfig.findGenes()

            findBestSignatureGenes(deseqDF, signatureGeneConfig)
                extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
        '''
        self.logger.info(f'BEGIN')
        
        bestDF = super().findGenes(deseqDF, fileName)

# example of bestDF
#          name     baseMean  log2FoldChange     lfcSE       stat         pvalue           padj
# 6       C1orf61  2958.098269       -7.635107  0.380337 -20.074565   1.231553e-89   3.179694e-86
# 3       ARPP21  1709.822182       -8.022019  0.399380 -20.086177   9.748529e-90   2.710541e-86
# 8       OLIG1  1323.934363       -7.620726  0.421037 -18.099881   3.193580e-73   2.355819e-70

        # select rows that are not poorDiscriminators
        selectRows = ~bestDF.loc[:, "name"].isin(self.poorDiscriminators)
        retDF = bestDF.loc[selectRows, :]

        self.logger.info(f'END')
        return retDF
    
    ################################################################################
    def _findPoorDiscriminators(self):
        self.logger.info("BEGIN")

        poorDiscriminators = []
        for key, items in self.intersectionDict.items():
            degree = len(key)
            if degree > self.degreeThreshold:
                poorDiscriminators = poorDiscriminators + items

        # return a sorted list of unique elements
        ret = sorted( list(set(poorDiscriminators)) )
        self.logger.info(f'len(poorDiscriminators) {len(poorDiscriminators)} poorDiscriminators :\n{ret}')

        self.logger.info("END")
        return ret

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

    cli = BestRemoveHighDegreeSignatureGeneConfigCLI(version=__version__ , 
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
    intersectionDictPath = cli.args.intersectionDict
    degreeThreshold = cli.args.degreeThreshold

    logger.warning(f'END')

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
                                 "--localCacheRoot", "myCacheRoot",
                                 "--intersectionDict", "./analysis/test/data/intersection.dict",
                                 "--degreeThreshold", "5"
                                ])


