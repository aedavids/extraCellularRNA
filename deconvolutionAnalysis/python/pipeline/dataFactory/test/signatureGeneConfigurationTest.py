#
# SignatureGeneConfigurationTest.py
# unit test for SignatureGenes.py 
# and pipeline.dataFactory.utilities.runSelectGenesOfInterest
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

import logging
import pandas as pd
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

##############################################################################
class SignatureGeneConfigurationTest(SignatureGeneConfiguration):
    logger = logging.getLogger(__name__)

    ################################################################################
    def __init__(self, msg : str, 
                 localCacheDir : str, 
                 title : str="plot title TODO"):
        super().__init__(
            dataSetName="testDataSet",
            design="~gender + category",
            padjThreshold=0.01,
            lfcThreshold=2.0,
            n=3,
            localCacheRootPath=localCacheDir,
            title=title
            )
        self.logger.info(f'BEGIN')
        self.msg = msg
        self.logger.info(f'END')

    ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        self.logger.info(f'BEGIN')
        self.logger.info(f'padjThreshold : {self.padjThreshold}')
        self.logger.info(f'msg : {self.msg}' )
        self.logger.info(f'n : {self.n}' )
        retDF = deseqDF.iloc[0:self.n, :]
        self.logger.info(f'END')
        return retDF
    
    ################################################################################
    def testDynamicalyloadedModules(self) -> str:
        return "hello world"
