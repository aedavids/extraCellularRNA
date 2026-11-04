#
# exampleCreateSignatureGeneConfig.py
# demonstrates how to implement the createSignatureGeneConfig() required
# by pipeline.py
# 
# Andrew E. Davidson
# aedavids@ucsc.edu
# 11/13/23
#

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2024-03-22'
__updated__ = '2024-03-22'

import logging

from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration
# from pipeline.dataFactory.test.signatureGeneConfigurationTest import SignatureGeneConfigurationTest
from analysis.bestLFCSignatureGeneConfig import BestLFCSignatureGeneConfig
from analysis.bestSignatureGeneConfigCLI import BestSignatureGeneConfigCommandLine

################################################################################
def createSignatureGeneConfig(outDirPath : str, vargs : list = None ) -> SignatureGeneConfiguration:
    '''
    
    arguments:
        outDirPath : str

        vargs : list
            option list of user defined command line arguments

    returns object derived from SignatureGeneConfiguration
    '''
    logger = logging.getLogger(__name__)

    logger.info("BEGIN")
    logger.info(f'outDirPath : {outDirPath}')

    if vargs :
        logger.info(f'optional CLI arguments vargs : {vargs}')
        tmsg = vargs[1]
    else :
        logger.info('no option command line arguments')
        tmsg = "integration unit test"


    cli = BestSignatureGeneConfigCommandLine(version=__version__ , 
                            author=__author__ ,
                            date=__date__, 
                            update=__updated__)
    
    cli.parse( vargs )

    logger.warning(f'command line arguments : {cli.args}')

    design          = cli.args.design 
    padjThreshold   = cli.args.padjThreshold
    lfcThreshold    = cli.args.lfcThreshold
    dataSetName     = cli.args.dataSetName
    number          = cli.args.number
    title           = cli.args.title
    localCacheRoot  = cli.args.localCacheRoot

    ret  = BestLFCSignatureGeneConfig(
                    dataSetName=dataSetName, 
                    design=design, 
                    padjThreshold=padjThreshold, 
                    lfcThreshold=lfcThreshold,
                    n=number, 
                    localCacheRootPath=localCacheRoot, 
                    title=title
    )

    logger.info("END")
    return ret
