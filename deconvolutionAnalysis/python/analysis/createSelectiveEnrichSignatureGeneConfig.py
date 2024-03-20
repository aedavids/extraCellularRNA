#
# createSelectiveEnrichSignatureGeneConfig.py
# ref: exampleCreateSignatureGeneConfig.py
# demonstrates how to implement the createSignatureGeneConfig() required
# by pipeline.py
# 
# Andrew E. Davidson
# aedavids@ucsc.edu
# 12/27/23
#

import logging

from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

from analysis.selectiveEnrichSignatureGeneConfig import SelectiveEnrichSignatureGeneConfig
from analysis.selectiveEnrichSignatureGeneConfigCLI import SelectiveEnrichSignatureGeneConfigCLI 
from analysis.selectiveEnrichSignatureGeneConfig import  __author__, __date__, __version__, __updated__

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


    cli = SelectiveEnrichSignatureGeneConfigCLI(version=__version__ , 
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
    intersectionDictPath = cli.args.intersectionDict
    numAdd          = cli.args.numAdd
    categories      = cli.args.classes


    ret  = SelectiveEnrichSignatureGeneConfig(
                    dataSetName=dataSetName, 
                    design=design, 
                    padjThreshold=padjThreshold, 
                    lfcThreshold=lfcThreshold,
                    n=number, 
                    localCacheRootPath=localCacheRoot, 
                    title=title,
                    intersectionDictionaryPath=intersectionDictPath,
                    numberOfGenesToAdd=numAdd,
                    categories=categories,
    )

    logger.info("END")
    return ret
