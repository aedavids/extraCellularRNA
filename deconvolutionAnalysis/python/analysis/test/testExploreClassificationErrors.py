#
# testExploreClassificationErrors.py
# unit test for exploreClassificationErrors.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/fractionsAsMulticlassClassification.ipynb
#


# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging

import numpy as np
import os
import pandas as pd
import pathlib as pl
import pprint as pp
import unittest

from analysis.exploreClassificationErrors import ExploreClassificationErrors
from analysis.metrics import CibersortResultsAsKWayClassifier
import analysis.test.testUtilities as tu

# AEDWIP_ENV = os.environ['AEDWIP_ENV']
# print(f'AEDWIP AEDWIP_ENV : {AEDWIP_ENV}')

ORIG_PYTHONPATH = os.environ['PYTHONPATH']
print(f'PYTHONPATH : {ORIG_PYTHONPATH}')
print(f'pwd: {os.getcwd()}')



################################################################################
class testExploreClassificationErrors(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger(__name__)

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
    def testFindFalseNegative(self):
        self.logger.info("BEGIN")

        evaluator = tu.getEvaluator(self.relativeRootPath)
        exploreErrors = ExploreClassificationErrors(evaluator)
        category = "BBB"
        falseNegativeSeries, groupByCountsSeries = exploreErrors.findFalseNegative(category)
        self.logger.info(f'category: {category} false negatives: \n{falseNegativeSeries}')
        self.logger.info(f'category: {category} groupByCountsSeries: \n{groupByCountsSeries}')

        expectedFalseNegativeSeries = pd.Series(["AAA", "FFF"], index=["s1", "s8"], name="predictedCat")
        expectedFalseNegativeSeries.index.name="Mixture"
        pd.testing.assert_series_equal(falseNegativeSeries, expectedFalseNegativeSeries)

        expectedGroupByCountsSeries = pd.Series([1,1], index=["AAA", "FFF"], name="predictedCat")
        expectedGroupByCountsSeries.index.name="predictedCat"
        pd.testing.assert_series_equal(groupByCountsSeries, expectedGroupByCountsSeries)

        self.logger.info("END\n")            

   ################################################################################
    def testFindFalsePositive(self):
        self.logger.info("BEGIN")

        evaluator = tu.getEvaluator(self.relativeRootPath)
        exploreErrors = ExploreClassificationErrors(evaluator)
        category = "BBB"
        falsePosSeries, groupByCountsSeries = exploreErrors.findFalsePositives(category)
        self.logger.info(f'category: {category} falsePosSeries: \n{falsePosSeries}')
        self.logger.info(f'category: {category} groupByCountsSeries: \n{groupByCountsSeries}')

        expectedFalsePosSeries = pd.Series(["EEE"], index=["s7"], name="category")
        expectedFalsePosSeries.index.name = 'Mixture'
        pd.testing.assert_series_equal(falsePosSeries, expectedFalsePosSeries)

        expectedGroupByCountsSeries = pd.Series([1], index=["EEE"], name="category")
        expectedGroupByCountsSeries.index.name="category"
        pd.testing.assert_series_equal(groupByCountsSeries, expectedGroupByCountsSeries)

        self.logger.info("END\n")            
                 
    ################################################################################
    def testFindSharedItems(self):
        self.logger.info("BEGIN")

        evaluator = tu.getEvaluator(self.relativeRootPath)
        exploreErrors = ExploreClassificationErrors(evaluator)

        intersectionDictPath = self.relativeRootPath.joinpath("data/intersection.dict")
        # 'Vagina_XXX_Whole_Blood': ['UVM_V_W'], 'UVM': ['UVM_AAA', 'UVM_V_W.1']
        setNames = set(['Vagina', 'Whole_Blood'])
        # find elements from all intersections containing Vagina and Whole_Blood
        intersectionsOfInterest = exploreErrors.findSharedGenes(intersectionDictPath, setNames)
        self.logger.info(f'intersectionsOfInterest:\n{intersectionsOfInterest}')
        expectedDict = {('Vagina', 'Whole_Blood'): ['UVM_V_W']}
        self.assertDictEqual(intersectionsOfInterest, expectedDict)

        # test degree = 1
        # find all intersections that contain elements from Whole_Blood
        setNames = set(['Whole_Blood'])
        intersectionsOfInterest = exploreErrors.findSharedGenes(intersectionDictPath, setNames)
        self.logger.info(f'degree = 1 intersectionsOfInterest:\n{intersectionsOfInterest}')

        expectedDict = {('Whole_Blood',): ['W_CCC', 'V_W'], 
                        ('Vagina', 'Whole_Blood'): ['UVM_V_W']}
        self.assertDictEqual(intersectionsOfInterest, expectedDict)        

        self.logger.info("END")

    
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
