#
# findMisclassificationErrors.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

# findMisclassificationErrors CommandLine display the doc string 
'''
    saves a data frame in csv format with  columns 'trueCat', 'predCat', 'errorCount'
'''

from argparse import ArgumentParser
#https://stackoverflow.com/a/3853776/4586180
from argparse import RawDescriptionHelpFormatter
from kimLabUtils.baseCommandLine import BBaseCommandLine

import logging
import pandas as pd
import os

# global variables use by FractionsCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-11-13'
__updated__ = '2023-11-13'

###############################################################################
class FindMisclassificationErrorsCLI( BBaseCommandLine ):
    '''
    Handle the command line, usage and help requests.
    '''
    ###############################################################################
    def __init__( self, version, author, date, update ):
        '''
        Implement a parser to interpret the command line argv string using argparse.

        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''
        super().__init__(version, author, date, update )

    ###############################################################################
    def _build( self ):
        self.parser = ArgumentParser( description=self._getLicence(), formatter_class=RawDescriptionHelpFormatter )
    
        #
        # optional arguments
        #

        # self.parser.add_argument( '-p', '--prefix', default=None, metavar="", type=str, action="store", 
        #                           help="save output files with prefix")

        self.parser.add_argument( '-t', '--threshold', default="1", metavar="",
                                            action='store', 
                                            type=int,
                                            help="select error counts >= threshold"
                                            + "\n default == 1"
        )     

        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )
        
        self.requiredArg.add_argument( '-c', '--confusionMatrixPath', required=True, default=".", metavar="",
                                            action='store', 
                                            help="path to csv file "
        ) 


        self.requiredArg.add_argument( '-o', '--outDir', required=True, default=".", metavar="",
                                            action='store', 
                                            help="locaiton to write output files"
        )
 
###############################################################################
def findMisclassificationErrors(
            cmDF : pd.DataFrame,
            threshold : int,
            ) -> pd.DataFrame:
    '''
    returns a data frame with three columns 'trueCat', 'predCat', 'errorCount'

    arguments
        cmDF : confusion matrix
            
        threshold : int
            find errors >= threshold
    '''
    
    # Find misclassifications greater than 'n'
    selectCols = ~cmDF.columns.isin(['true\predicted'])
    cmValuesDF = cmDF.loc[:, selectCols]
    misclassificationsLogicalDF = (cmValuesDF >= threshold) & (cmValuesDF != cmValuesDF.values.diagonal())
    # #print(f'misclassifications \n{misclassifications}')
    
    # # Get row and column indices
    row_indices, col_indices = misclassificationsLogicalDF.to_numpy().nonzero()
    
    trueCategories = cmDF.loc[:, "true\predicted"].tolist()
    predictedCategories = cmDF.columns[selectCols].tolist()
    
    rowNames = []
    colNames = []
    errorCounts = []
    for row, col in zip(row_indices, col_indices):
        #print(f"{trueCategories[row]}, {predictedCategories[col]}, Error Count: {cmValuesDF.iloc[row, col]}")    
        rowNames.append( trueCategories[row] )
        colNames.append( predictedCategories[col] )
        errorCounts.append( cmValuesDF.iloc[row, col] )
    
    misClassificationDF = pd.DataFrame({ 
                                        "trueCat" : rowNames,
                                        "predCat" : colNames,
                                        "errorCount" : errorCounts
                                })    

    return misClassificationDF.sort_values(by="errorCount", ascending=False)

################################################################################
def main(inCommandLineArgsList=None):
    '''
    ref: pipeline/dataFactory/test/testInjection.py 

    ```
    see extraCellularRNA/deconvolutionAnalysis/bin/pipeline.sh
    ```

    see unit test unit test 
    extraCellularRNA/deconvolutionAnalysis/python/analysis/test/testFractions.py
    
    '''
    # we only configure logging in main module
    loglevel = "WARN"
    #loglevel = "INFO"
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    logger = logging.getLogger(__file__)
    logger.warning("BEGIN")

    #
    # always log run time env to make debug easier
    #
    ORIG_PYTHONPATH = os.environ['PYTHONPATH']
    logger.warning(f'PYTHONPATH : {ORIG_PYTHONPATH}')
    logger.warning(f'FILE: {__file__}')
    logger.warning(f'PWD: {os.getcwd()}')

    #
    # parse command line arguments
    #
    cli = FindMisclassificationErrorsCLI( 
                    version=__version__ , 
                    author=__author__ ,
                    date=__date__, 
                    update=__updated__)

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    confusionMatrixPath = cli.args.confusionMatrixPath
    threshold           = cli.args.threshold
    outDir              = cli.args.outDir

    logger.warning(f' confusionMatrixPath: {confusionMatrixPath}')
    logger.warning(f' threshold : {threshold}')
    logger.warning(f' outDir : {outDir}')

    os.makedirs(outDir, exist_ok=True)

    cmDF = pd.read_csv(confusionMatrixPath)
    retDF = findMisclassificationErrors(cmDF, threshold)
    
    if outDir[-1] == "/":
        outDir = outDir[0:-1]

    savePath = f'{outDir}/classificationErrors.csv'
    retDF.to_csv(savePath)
    logger.warning(f'saved classificationErrors.csv to : {savePath}')

################################################################################
if __name__ == '__main__':
    main()

    # main( inCommandLineArgsList=["--confusionMatrixPath", "my_confusionMatrixPath",
    #                                 "--outDir", "my_outDir",
    #    ])
    
    #main(["-h"])