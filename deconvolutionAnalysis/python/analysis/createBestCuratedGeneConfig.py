#
# createBestCuratedGeneConfig.py
# ref: exampleCreateSignatureGeneConfig.py
# demonstrates how to implement the createSignatureGeneConfig() required
# by pipeline.py
# 
# Andrew E. Davidson
# aedavids@ucsc.edu
# 01/08/24
#

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-11-13'
__updated__ = '2023-11-13'

import logging

from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration
from analysis.bestCuratedGeneConfig import BestCuratedGeneConfig
from analysis.bestCuratedGeneConfigCLI import BestCuratedGeneConfigCommandLine

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


    cli = BestCuratedGeneConfigCommandLine(version=__version__ , 
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
    interesectionDictPath = cli.args.interesectionDictPath

    if not cli.args.ascending:
        # backwards compatiblity default
        ascending = True
    else :
        ascending = False
    logger.info(f' ascending : {ascending}')

    
    ret  = BestCuratedGeneConfig(
                    dataSetName=dataSetName, 
                    design=design, 
                    padjThreshold=padjThreshold, 
                    lfcThreshold=lfcThreshold,
                    n=number, 
                    localCacheRootPath=localCacheRoot, 
                    title=title,
                    interesectionDictPath=interesectionDictPath,
                    ascending=ascending
    )

    logger.info("END")
    return ret
