#
# exampleCreateSignatureGeneConfig.py
# demonstrates how to implement the createSignatureGeneConfig() required
# by pipeline.py
# 
# Andrew E. Davidson
# aedavids@ucsc.edu
# 11/13/23
#

import logging

from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration
from pipeline.dataFactory.test.signatureGeneConfigurationTest import SignatureGeneConfigurationTest

def createSignatureGeneConfig(outDirPath : str, vargs : list = None ) -> SignatureGeneConfiguration:
    '''
    example of how to write a module that can be "dynamically injected" into our pipeline
    
    arguments:
        outDirPath : str

        vargs : list
            option list of user defined command line arguments

    returns object derived from SignatureGeneConfiguration
    '''
    logger = logging.getLogger("__name__")

    logger.info("BEGIN")
    logger.info(f'outDirPath : {outDirPath}')

    if vargs :
        logger.info(f'optional CLI arguments vargs : {vargs}')
        tmsg = vargs[1]
    else :
        logger.info('no option command line arguments')
        tmsg = "integration unit test"

    logger.info(f'tmsg : {tmsg}')
    ret  = SignatureGeneConfigurationTest(msg=tmsg, localCacheDir=outDirPath )

    logger.info("END")
    return ret
