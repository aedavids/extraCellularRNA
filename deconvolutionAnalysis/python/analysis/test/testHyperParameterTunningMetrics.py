#
# testHyperParameterTunningMetrics.py
# unit test for HyperParameterTunningMetrics.py
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

# from analysis.metrics import CibersortResultsAsKWayClassifier
# import analysis.test.testUtilities as tu
import numpy as np
import pandas as pd
import unittest

from analysis.hyperParameterTunningMetrics import adjacentRowSort

################################################################################
class TestFractions(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger(os.path.basename(__file__))

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
    def testAdjacentRowSort(self):
        '''
        make sure we can load data correctly
        '''
        self.logger.info("BEGIN")

        inputDict = {
            "trueCat" :["Skin_Sun_Exposed_Lower_leg", "Colon_Transverse", "Breast_Mammary_Tissue",
                        "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Muscularis", 
                        "Esophagus_Gastroesophageal_Junction"], 

            "predCat" :["Skin_Not_Sun_Exposed_Suprapubic", "Colon_Sigmoid", "Adipose_Subcutaneous", 
                        "Skin_Sun_Exposed_Lower_leg", "Esophagus_Gastroesophageal_Junction", 
                        "Esophagus_Muscularis"],

            "numericCol" : [116, 90, 84, 82, 67, 66]
        }
        
        df = pd.DataFrame( inputDict )

        retDF = adjacentRowSort(df, "trueCat", "predCat", "numericCol", verbose=False)
        self.logger.info(f'retDF:\n{retDF}')

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
