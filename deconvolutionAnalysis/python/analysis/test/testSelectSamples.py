#
# testSelectSamples.py.py
# unit test for selectSamples.py.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/fractionsAsMulticlassClassification.ipynb
#


# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging

# import numpy as np
import os
import pandas as pd
import pathlib as pl
# import pprint as pp
import unittest


from analysis.selectSamples import SelectSamples

# AEDWIP_ENV = os.environ['AEDWIP_ENV']
# print(f'AEDWIP AEDWIP_ENV : {AEDWIP_ENV}')

ORIG_PYTHONPATH = os.environ['PYTHONPATH']
print(f'PYTHONPATH : {ORIG_PYTHONPATH}')
print(f'pwd: {os.getcwd()}')

################################################################################
class testSelectSamples(unittest.TestCase):
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
    def testSelect(self):
        self.logger.info("BEGIN")

        resultsFilePath = self.relativeRootPath.joinpath("data/results.tsv")
        ss = SelectSamples(resultsFilePath, separator='\t')

        samples = ['s1', 's3', 's5']
        samplesDF = ss.select( samples )
        self.logger.info(f'retDF:\n{samplesDF}')


        expectedPath = self.relativeRootPath.joinpath("data/testSelectSamples/selectSamples.csv")
        expectedDF = pd.read_csv(expectedPath, index_col="Mixture")
        self.logger.info(f'expectedDF:\n{expectedDF}')

        pd.testing.assert_frame_equal(expectedDF, samplesDF)

        self.logger.info("END\n")            

   ################################################################################
    def testSimplify(self):
        self.logger.info("BEGIN")

        resultsFilePath = self.relativeRootPath.joinpath("data/results.tsv")
        ss = SelectSamples(resultsFilePath, separator='\t')

        samples = ['s1', 's3', 's5']
        samplesDF = ss.select( samples )
        self.logger.info(f'samplesDF:\n{samplesDF}')

        simplifiedDF = ss.simplify(samplesDF, threshold=0.3)
        self.logger.info(f'simplifiedDF:\n{simplifiedDF}')

        expectedPath = self.relativeRootPath.joinpath("data/testSelectSamples/simplifiedSamples.csv")
        #simplifiedDF.to_csv(expectedPath)
        expectedDF = pd.read_csv(expectedPath,index_col="Mixture")

        pd.testing.assert_frame_equal(expectedDF, simplifiedDF)

        self.logger.info("END\n")            
  

    
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
