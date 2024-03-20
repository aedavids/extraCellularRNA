#
# exploreClassificationErrors.py
# functions for analysising the results from CIBERSORTx as if it was a multi-label classifier
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ExploreClassificationErrorsCLI CommandLine display the doc string 
'''
    TODO doc string
    See analysis.test.testExploreClassificationErrors
'''

import ast 
import logging
import numpy as np
import os
import pandas as pd
import pprint as pp

from analysis.metrics import CibersortResultsAsKWayClassifier
from analysis.exploreClassificationErrorsCLI import ExploreClassificationErrorsCLI

# global variables use by ExploreClassificationErrorsCLI
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-11-13'
__updated__ = '2023-11-13'

##############################################################################
class ExploreClassificationErrors( object ):
    '''
    TODO

    main() # cli
    findFalseNegative()
    findSharedGenes()
    findFalsePositives()
  
    '''

    logger = logging.getLogger(__name__)

    ###############################################################################
    def __init__( self,  
                    evaluator : CibersortResultsAsKWayClassifier):
        '''
        TODO
        '''
        self.evaluator = evaluator

    ###############################################################################
    def findFalseNegative( self, category : str )  -> tuple[pd.Series, pd.Series] :
        '''
        ref: analysis.testExploreClassificationErrors testFindFalseNegative()

        arguments:
            category : str

        returns:
            (falseNegativeSeries, groupByCountsSeries)

            falseNegativeSeries : 
                index is the example id. The values is the predicted category

            groupByCountsSeries :
                index is the predicted category
                the values is the number of examples that where incorreclty predicted
        '''

        # get all the sample that are of type category
        expectedExamplesSeriesRows = self.evaluator.trueLablesSeries == category

        # get the predictions for these samples
        pDF = self.evaluator.fractionsWithPredictionsDF
        predictitionsSeries = pDF.loc[expectedExamplesSeriesRows,'predictedCat']

        # find the false negative
        falseNegativeSeriesRows = predictitionsSeries != category
        falseNegativeSeries = predictitionsSeries.loc[falseNegativeSeriesRows]

        # group by the values in the series
        groupByCountsSeries = falseNegativeSeries.groupby( falseNegativeSeries ).count()

        return (falseNegativeSeries, groupByCountsSeries)

    ###############################################################################
    def findFalsePositives( self, category : str) -> tuple[pd.Series, pd.Series] :
        '''
        ref: analysis.testExploreClassificationErrors testFindFalsePositive()
        '''

        # get all examples that where prediced to be of type category
        pDF = self.evaluator.fractionsWithPredictionsDF
        predictitionsSeriesRows = pDF.loc[:,"predictedCat"] == category
        predictitionsSeries = pDF.loc[predictitionsSeriesRows, 'predictedCat']

        # find truth value for examples that where predicted to be of type category
        truthSeries  = self.evaluator.trueLablesSeries.loc[predictitionsSeries.index]

        # compare predicted to truth
        falsePosRows = predictitionsSeries != truthSeries
        falsePosSeries = truthSeries.loc[falsePosRows]
        groupByCountsSeries = falsePosSeries.groupby( falsePosSeries ).count()

        return (falsePosSeries, groupByCountsSeries)
      
    ###############################################################################
    def findSharedGenes(self, intersectionDictPath : str, 
                            setNames : set[str]
                            ) -> dict :
        '''
        TODO

        ref: analysis.testExploreClassificationErrors testFindSharedItems()
        '''
        retDict = dict()
        setNamesSet = set(setNames)

        with open(intersectionDictPath) as f: 
            data = f.read() 

        # we can not convert the intersectionData to a DataFrame
        # the size of the interesections is not constant. 
        # we would need to pad the dictionary values
        intersectionDict = ast.literal_eval(data)

        # find intersections that contain sets of interested
        # we do not know what order the set names are in  
        degree = len(setNames)
        for keyTuple, values in intersectionDict.items():
            #self.logger.error(f'AEDWIP debug keyTuple: {keyTuple}')

            #hack to make debug eaiser
            # if keyTuple[0] == 'ACC' and keyTuple[1] ==  'Adipose_Subcutaneous' and keyTuple[2] ==  'Adipose_Visceral_Omentum':
            #     self.logger.error(f'AEDWIP hack break point')

            if len(keyTuple) == 0:
                # empty set is always a subset of all other sets
                continue

            if len(keyTuple) < degree:
                # degree is the number of sets that have elements in the intersection
                # we want all interesection that contain elements from setNames
                continue

            s = set(keyTuple)
            # bug if s.issubset( setNamesSet ) :
            if setNamesSet.issubset( s ) :
                saveKey = tuple( sorted(s) )
                retDict[saveKey] = values

        return retDict

################################################################################
def main(inCommandLineArgsList=None):
    '''
    
    ```

    see unit test unit test 
    analysis.test.testExploreClassificationErrors
    
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
    cli = ExploreClassificationErrorsCLI( 
                    version=__version__ , 
                    author=__author__ ,
                    date=__date__, 
                    update=__updated__)


    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    #logger.error(f'AEDWIP cli: {cli}')
    
    cmd = cli.args.subcommand
    logger.debug(f'AEDWIP cmd: {cmd}')

    outDir                = cli.args.outDir
    os.makedirs(outDir, exist_ok=True)


    if cmd == "fn" or cmd == "fp":
        category              = cli.args.category
        resultsPath           = cli.args.results
        expectedFractionsPath = cli.args.expected
        evaluator = CibersortResultsAsKWayClassifier(resultsPath, expectedFractionsPath)
        ece = ExploreClassificationErrors(evaluator)

        if cmd == "fn":
            logger.debug(f'AEDWIP if cmd fn')
            falseSeries, groupByCountsSeries = ece.findFalseNegative(category)
            fileName = f"{outDir}/{category}.falseNegative.csv"
            groupByFileName=f"{outDir}/{category}.falseNegativeGroupByCounts.csv"


        else:
            logger.debug(f'AEDWIP if cmd fb')
            falseSeries, groupByCountsSeries = ece.findFalsePositives(category)
            fileName = f"{outDir}/{category}.falsePositive.csv"
            groupByFileName=f"{outDir}/{category}.falsePositiveGroupByCounts.csv"

        falseSeries.to_csv(fileName)
        logger.warning(f'saving file: {fileName}')
        logger.warning(f'saving file: {groupByFileName}')
        groupByCountsSeries.to_csv(groupByFileName)

    else :
        # cmd = sq
        ece = ExploreClassificationErrors(None) # hack
        intersectionDictFile = cli.args.intersectionDictionary

        categories           = cli.args.category
        retDict = ece.findSharedGenes(intersectionDictFile, categories)

        setNames = "-".join(categories)
        fileName = f"{outDir}/{setNames}.sharedItems.dict"
        logger.warning(f'saving file: {fileName}')
        with open (fileName, "w") as fd:
            fd.write( pp.pformat(retDict, indent=4, sort_dicts=True))


################################################################################
if __name__ == '__main__':
    main()
