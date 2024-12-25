"""
testCreateCiberSortMixtureMatrix.py

Unit tests for the reateCiberSortMixtureMatrix module.

Author: Andrew E. Davidson
aedavids@ucsc.edu

Created: 20224-12-25
"""

import os
import pandas as pd
import pprint as pp
# AEDWIP_ENV = os.environ['AEDWIP_ENV']
# print(f'AEDWIP AEDWIP_ENV : {AEDWIP_ENV}')

ORIG_PYTHONPATH = os.environ['PYTHONPATH']
print(f'PYTHONPATH : {ORIG_PYTHONPATH}')
print(f'FILE: {__file__}')
print(f'PWD: {os.getcwd()}')

# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging

import numpy as np
import pandas as pd
import pathlib as pl
import unittest

from tempus.createCiberSortMixtureMatrix import createMixtureMatrix

################################################################################
class TestCreateCiberSortMixtureMatrix(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger(__name__)


    ################################################################################
    def _getNormalizedCountData(self) -> pd.DataFrame :
        '''
        Get the tempus, "complete seq" test count data used to create the signature matrix
        '''
        self.logger.info(f'BEGIN')
       
        data = [
                ['(A)n',       '(A)n',       'Microsatellite', 1, 2, 3],
                ['(AAA)n',     '(AAA)n',     'Microsatellite', 4, 5, 6],
                ['(AAAAAAC)n', '(AAAAAAC)n', 'Microsatellite', 7, 8, 9],
                ['(AAAAAAG)n', '(AAAAAAG)n', 'Microsatellite',  10, 11, 12],
                ['(AAAAAAT)n', '(AAAAAAT)n', 'Microsatellite', 13, 14, 15]
                ]

        cols = ['gene_id', 'gene_name', 'gene_biotype', 'SLDK3_T1_100_S1_L007', 'SLDK3_T1_100K_S2_L007', 'SLDK3_T1_10K_S3_L007']

        df = pd.DataFrame(data, columns=cols)
        self.logger.info(f'df\n{df}')
        self.logger.info(f'END')
        return df

    ################################################################################
    def testCreateMixtureMatrix(self):
        '''
        '''
        self.logger.info(f'BEGIN')

        countDF = self._getNormalizedCountData()
        countDF.set_index('gene_id', inplace=True)

        # drop the gene_name and biotype columns
        countDF = countDF.drop(columns=["gene_name", "gene_biotype"])

        self.logger.info(f'countData\n{countDF}')

        #createMixtureMatrix(countDF, genesOfInterest = ["X7D_LINE", "Zaphod"])

        self.logger.info(f'END')

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
