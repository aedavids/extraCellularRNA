#
# createHistoryGeneSet.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
#

import logging
import os
import pandas as pd
import pprint as pp
import sys

from analysis.utilities import findFile
from pipeline.dataFactory.driver import _countExtraHeaderLines

################################################################################
def createHistoryGeneSet(resultDirs : list[str]):
    '''
    select all *"*vs_all.results" in hte list of resultDirs
    and creates a geneSet from the 'name' column

    arguments:
        resultsDirs
            a list of directory paths that contain deseq results
    
    print geneSet to stdout
    '''
    pass

################################################################################
def main():
    '''
    '''
    # we only configure logging in main module
    loglevel = "WARN"
    #loglevel = "INFO"
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

    logger.error(f'Number of arguments: {len(sys.argv)}, arguments.')
    argv = sys.argv
    logger.error(f'Argument List:, {str(sys.argv)}' )
    logger.error(f'Argument List:, {sys.argv}' )


    if len(argv) <= 1:
        #print help
        print('\nERROR missing list of directories')
        print(createHistoryGeneSet.__doc__)
        sys.exit(1)

    retGeneSet = set()
    for i in range(1, len(argv)):
        dirPath = argv[i]
        logger.info(f'i: {i} dirPath: {dirPath}')
        pattern="*vs_all.results"
        resultsFile = findFile(dirPath, pattern)
        for filePath in resultsFile:
            numRowsToSkip = _countExtraHeaderLines(filePath)            
            df = pd.read_csv(filePath, header=numRowsToSkip)
            genes = df.loc[:, 'name'].to_list()
            retGeneSet = retGeneSet.union( set(genes) )

    print( pp.pformat(retGeneSet) )
    logger.info(f'where does log info go?')
    #print('print goest to stdout log does not')

    logger.warning(f'END')

################################################################################
if __name__ == '__main__':
    '''
    hack used to test parsing of vargs
    '''
    # /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/analysis/test/data/testEnrichSignatureGeneConfig/1vsAll
    # /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/analysis/test/data/testSelectiveEnrichSignatureGeneConfig/1vsAll
    main()
