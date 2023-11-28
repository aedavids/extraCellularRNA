#
# testFractions.py
# unit test for fractions.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/fractionsAsMulticlassClassification.ipynb
#

import os
import pandas as pd
import pathlib as pl
import pprint as pp
# AEDWIP_ENV = os.environ['AEDWIP_ENV']
# print(f'AEDWIP AEDWIP_ENV : {AEDWIP_ENV}')

ORIG_PYTHONPATH = os.environ['PYTHONPATH']
print(f'PYTHONPATH : {ORIG_PYTHONPATH}')
print(f'pwd: {os.getcwd()}')

# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging

from analysis.metrics import CibersortResultsAsKWayClassifier
import numpy as np
import pandas as pd
import unittest

################################################################################
class TestFractions(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger(__name__)

    # https://stackoverflow.com/a/14493895/4586180
    maxDiff = None

    #pwd = os.getcwd()

    # python -m unittest discover .
    # we need to be able to find test file relative to module location
    # os.getcwd() is often the top of our source code tree
    # https://stackoverflow.com/a/57836145/4586180
    relativeRootPath = pl.Path(os.path.dirname(os.path.realpath(__file__)))

    ################################################################################
    @unittest.skip("skip template. it is not a real test")
    def testTemplate(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n")            
                
    ################################################################################
    def testAAALoad(self):
        '''
        make sure we can load data correctly
        '''
        self.logger.info("BEGIN")

        # # dataDir = "../python/analysis/test/data"
        # dataDir = "./data"
        # expectedFractionsPath = f'{dataDir}/expectedFractions.tsv'
        # resultPath            = f'{dataDir}/results.tsv'

        # self.logger.info(f'pwd: {os.getcwd()}')
        # self.logger.info(f"expectedFractionsPath :{expectedFractionsPath}")

        # # evaluator is a object that calculates classification statistics
        # evaluator = CibersortResultsAsKWayClassifier(resultPath, expectedFractionsPath, verbose=True)

        evaluator = self._getEvaluator()
        # sanity test make sure our test data has not accidently changed
        trueLabelSeries = evaluator.trueLablesSeries
        self.logger.info(f'trueLablesSeries\n{trueLabelSeries}')
        
        predictedLabelsSeries = evaluator.predictedLabelsSeries
        self.logger.info(f'predictedLabelsSeries\n{predictedLabelsSeries}')

        expectedTrueLabelsSeries, expectedPredictedLabelsSeries = self._getExpectedTrueLabels(evaluator)
        
        pd.testing.assert_series_equal(expectedTrueLabelsSeries, trueLabelSeries)
        pd.testing.assert_series_equal(expectedPredictedLabelsSeries, predictedLabelsSeries)

        self.logger.info("END\n")            

    ################################################################################
    def testBBBConfusionMatrix(self):
        '''
        make sure the stats and error metrics are calculated as expected
        '''
        self.logger.info("BEGIN")

        # Note we have 4 classes AAA, BBB, CCC,DDD
        # There are no DDD examples or predictions. 
        # EEE has one expected, prediction is wrong. if we do not make a prediction join will drop EEE
        # FFF has no expected examples but one predicted
        evaluator = self._getEvaluator()
        #self.logger.info(f'evaluator.confusionDF.to_dict()\n{evaluator.confusionDF.to_dict()}')

        # 
        # test confussion matrix is calculated correctly for all possible class types
        #
        # true\predicted,AAA,BBB,CCC,DDD,EEE,FFF
        # AAA [[1 0 0 0 0 0] s1 true = AAA, s1 predicted = AAA
        # BBB [1 2 0 0 0 1] s2 s4 s5, S8 true = BBB s2 predicted AAA, S4 predicted BBB or CCC S5 prediced BBB S8 predicted FFF
        # CCC [0 0 2 0 0 0] s3, s6 true = CCC, s3 s6 predicted CCC
        # DDD [0 0 0 0 0 0] DDD no true or predicted
        # EEE [0 1 0 0 0 0] s7 true = EEE s7 predicted BBB
        # FFF [0 0 0 0 0 0]]] = FFF no TRUE, S8 predicted FFF
        #
        expectedLabels = ['AAA', 'BBB', 'CCC', 'DDD', 'EEE', 'FFF']
        # cm, cmReportStr, cmReportDict = evaluator.calculate(expectedLabels)
        cm, cmReportDF = evaluator.calculate(expectedLabels)
        self.logger.info(f'confusion matrix labels = expectedLabels\n{cm}')
        self.assertEqual(expectedLabels, evaluator.labels)
        
        # self.logger.info(f'cmReportStr:\n{cmReportStr}')
        self.logger.info(f'cmReportDict:')
        print(cmReportDF)
        # pp.pprint(cmReportDF.to_dict(), sort_dicts=False)
        expectedCmReportDict = self._getExpectedReportDict()
        
        expectedCmReportDF = pd.DataFrame(expectedCmReportDict)
        self.logger.info(f'\n********** AEDWIP pandas')
        print(expectedCmReportDF)
        # print(cmReportDF.transpose(copy=True))
        pd.testing.assert_frame_equal(expectedCmReportDF, cmReportDF)

        expectedCM = self._getExpectedConfusionMatrix()
        self.assertTrue( (expectedCM == cm).all() )

        #
        # test confusion matrix is calculated correctly when we only
        # want classes that that appear at least once in y_true or y_pred
        # 
        cmLablesNone, cmLablesNoneReportDF = evaluator.calculate(labels=None)
        #aedwip cmLablesNone = evaluator.confusionMatrix
        self.logger.info(f'confusion matrix labels = None \n{cmLablesNone}')

        expectedCMLablesNone = self._getExpectedCMLablesNone()
        self.assertTrue( (expectedCMLablesNone == cmLablesNone).all() )
        self.logger.info(f'expectedCMLablesNone:\n{expectedCMLablesNone}')
        self.logger.info(f'cmLablesNoneReportDF:\n{cmLablesNoneReportDF}')

        self.logger.info("END\n")            
   
    ################################################################################
    def testCCCConfusionMatrixSubset(self):
        '''
        make sure the stats and error metrics are calculated as expected
        '''
        self.logger.info("BEGIN")   
        evaluator = self._getEvaluator()
        #self.logger.info(f'evaluator.confusionDF.to_dict()\n{evaluator.confusionDF.to_dict()}')

        # 
        # test confussion matrix is calculated correctly for subset of classes class types
        #
        expectedLabels = ['AAA', 'BBB', 'CCC']
        # cm, cmReportStr, cmReportDict = evaluator.calculate(expectedLabels)
        cm, cmReportDF = evaluator.calculate(expectedLabels)
        
        # precision  recall  f1-score  support are the same as when we run with all possible labes
        # obviously tp, fn, fp, and tn change. we have few sample so counts are lower
        # specificity and sensitivity are not the same
        self.logger.info(f'cmReportDF:]\n{cmReportDF}')
        self.logger.info("END")   

    #
    # utility functions
    #

    ################################################################################
    def _getExpectedConfusionMatrix(self) :
        '''
                    predicted
                    A  B  C  D  E  F
            array([ A [1, 0, 0, 0, 0, 0],
                    B [1, 2, 0, 0, 0, 1],
                    C [0, 0, 2, 0, 0, 0],
            true    D [0, 0, 0, 0, 0, 0],
                    E [0, 1, 0, 0, 0, 0],
                    F [0, 0, 0, 0, 0, 0]])        
        '''
                                   # predicted
                                   # A  B  C  D  E  F
        expectedCM =  np.array( [   [1, 0, 0, 0, 0, 0],
                                    [1, 2, 0, 0, 0, 1],
                                    [0, 0, 2, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0],
                                    [0, 1, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0]] )
        
        return expectedCM
    
    ################################################################################
    def _getExpectedCMLablesNone(self):
        '''
        self.confusionMatrix
                predicted
                    A  B  C  E  F
        array([   A [1, 0, 0, 0, 0],
        true      B [1, 2, 0, 0, 1],
                C [0, 0, 2, 0, 0],
                E [0, 1, 0, 0, 0],
                F [0, 0, 0, 0, 0]])        
        '''

                                           # predicted
                                           # A  B  C  E  F
        expectedCMLablesNone = np.array( [ [ 1, 0, 0, 0, 0],
                                            [1, 2, 0, 0, 1],
                                            [0, 0, 2, 0, 0],
                                            [0, 1, 0, 0, 0],
                                            [0, 0, 0, 0, 0]] )       
        return expectedCMLablesNone 

    ################################################################################
    def _getEvaluator(self) -> CibersortResultsAsKWayClassifier:
        
        #dataDir = os.path.join(self.pwd, "../python/analysis/test/data")
        #dataDir = os.path.join(self.pwd, "./data")
        dataDir = self.relativeRootPath.joinpath("data")
        # expectedFractionsPath = f'{dataDir}/expectedFractions.tsv'
        # resultPath            = f'{dataDir}/results.tsv'
        expectedFractionsPath = dataDir.joinpath('expectedFractions.tsv')
        resultPath            = dataDir.joinpath('results.tsv')

        self.logger.info(f'pwd: {os.getcwd()}')
        self.logger.info(f"expectedFractionsPath :{expectedFractionsPath}")

        # evaluator is a object that calculates classification statistics
        evaluator = CibersortResultsAsKWayClassifier(resultPath, expectedFractionsPath, verbose=True)

        return evaluator

    ################################################################################
    def _getExpectedTrueLabels(self, 
                               evaluator : CibersortResultsAsKWayClassifier) -> tuple[pd.Series, pd.Series]:
        trueLabelDict = {
        's1' : 'BBB',
        's2' : 'AAA',
        's3' : 'CCC',
        's4' : 'BBB',
        's5':'BBB',
        's6': 'CCC',
        's7': 'EEE',
        's8': 'BBB'}
        trueLabelsSeries = pd.Series(trueLabelDict, name="category")
        trueLabelsSeries.index.name = "Mixture"   

        predictedLabeslDict = {
        's1' : 'AAA',
        's2' : 'AAA',
        's3' : 'CCC',
        's4' : 'BBB',
        's5' : 'BBB',
        's6' : 'CCC',
        's7' : 'BBB',
        's8' : 'FFF'}
        
        predictedLabelsSeries = pd.Series(predictedLabeslDict, name="predictedCat")
        predictedLabelsSeries.index.name = "Mixture"   

        ret = (trueLabelsSeries, predictedLabelsSeries)
        return ret

    ################################################################################
    def _getExpectedReportDict(self):
        retDict = {
            'precision': {'AAA': 0.5,
               'BBB': 0.6666666666666666,
               'CCC': 1.0,
               'DDD': np.nan,
               'EEE': np.nan,
               'FFF': 0.0,
               'micro avg': 0.625,
               'macro avg': 0.5416666666666666,
               'weighted avg': 0.738095238095238},
            'recall': {'AAA': 1.0,
                        'BBB': 0.5,
                        'CCC': 1.0,
                        'DDD': np.nan,
                        'EEE': 0.0,
                        'FFF': np.nan,
                        'micro avg': 0.625,
                        'macro avg': 0.625,
                        'weighted avg': 0.625},
            'f1-score': {'AAA': 0.6666666666666666,
                        'BBB': 0.5714285714285715,
                        'CCC': 1.0,
                        'DDD': np.nan,
                        'EEE': np.nan,
                        'FFF': np.nan,
                        'micro avg': 0.625,
                        'macro avg': 0.746031746031746,
                        'weighted avg': 0.7074829931972789},
            'support': {'AAA': 1.0,
                        'BBB': 4.0,
                        'CCC': 2.0,
                        'DDD': 0.0,
                        'EEE': 1.0,
                        'FFF': 0.0,
                        'micro avg': 8.0,
                        'macro avg': 8.0,
                        'weighted avg': 8.0},
            'specificity': {'AAA': 0.8571428571428571,
                            'BBB': 0.75,
                            'CCC': 1.0,
                            'DDD': 1.0,
                            'EEE': 1.0,
                            'FFF': 0.875,
                            'micro avg': np.nan,
                            'macro avg': np.nan,
                            'weighted avg': np.nan},
            'sensitivity': {'AAA': 1.0,
                            'BBB': 0.5,
                            'CCC': 1.0,
                            'DDD': np.nan,
                            'EEE': 0.0,
                            'FFF': np.nan,
                            'micro avg': np.nan,
                            'macro avg': np.nan,
                            'weighted avg': np.nan},
            'tp': {'AAA': 1.0,
                    'BBB': 2.0,
                    'CCC': 2.0,
                    'DDD': 0.0,
                    'EEE': 0.0,
                    'FFF': 0.0,
                    'micro avg': np.nan,
                    'macro avg': np.nan,
                    'weighted avg': np.nan},
            'fn': {'AAA': 0.0,
                    'BBB': 2.0,
                    'CCC': 0.0,
                    'DDD': 0.0,
                    'EEE': 1.0,
                    'FFF': 0.0,
                    'micro avg': np.nan,
                    'macro avg': np.nan,
                    'weighted avg': np.nan},
            'fp': {'AAA': 1.0,
                    'BBB': 1.0,
                    'CCC': 0.0,
                    'DDD': 0.0,
                    'EEE': 0.0,
                    'FFF': 1.0,
                    'micro avg': np.nan,
                    'macro avg': np.nan,
                    'weighted avg': np.nan},
            'tn': {'AAA': 6.0,
                    'BBB': 3.0,
                    'CCC': 6.0,
                    'DDD': 8.0,
                    'EEE': 7.0,
                    'FFF': 7.0,
                    'micro avg': np.nan,
                    'macro avg': np.nan,
                    'weighted avg': np.nan}}
        
        return retDict


################################################################################
if __name__ == "__main__":
    # we only configure logging in main module
    # loglevel = p.getProperty("LOG_LEVEL")
    loglevel = "INFO"
    # logFMT = p.getProperty("LOG_FMT")
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
