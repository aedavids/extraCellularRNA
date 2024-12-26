"""
testCreateCibersortSignature.py

Unit tests for the createCibersortSignatureMatrix module.

Author: Andrew E. Davidson
aedavids@ucsc.edu

Created: 20224-12-22
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

from tempus.createCibersortSignatureMatrix import createSignatureMatrix

################################################################################
class TestCreateCibersortSignatureMatrix(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger(__name__)


    # python -m unittest discover .
    # we need to be able to find test file relative to module location
    # os.getcwd() is often the top of our source code tree
    # https://stackoverflow.com/a/57836145/4586180
    relativeRootPath = pl.Path(os.path.dirname(os.path.realpath(__file__)))
    
    #expectedSignatureFilePath = "data/geneSignatureProfiles/best/ciberSortInput/signatureGenes.tsv"
    expectedSignatureFilePath = relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/signatureGenes.tsv")
    
    localCacheDir= relativeRootPath.joinpath("data/tmp")
    logger.info(f'localCacheDir: {localCacheDir}')

    ################################################################################
    def _getNormalizedCountData(self) -> pd.DataFrame :
        '''
        Get the tempus, "complete seq" test count data used to create the signature matrix
        '''
        self.logger.info(f'BEGIN')
        dict = {}
        dict[ "gene_id"               ] = ["X7D_LINE", "X8_LINE", "X9_LINE", "Zaphod", "Zaphod2", "Zaphod3" ]
        dict[ "gene_name"             ] = ["X7D_LINE",  "X8_LINE", "X9_LINE", "Zaphod", "Zaphod2", "Zaphod3" ]
        dict[ "biotype"               ] = ["LINE",     "LINE",     "LINE",    "DNA",    "DNA",     "DNA"]
        dict[ "SLDK3_T1_100_S1_L007"  ] = [ 1.0,  2.0,  3.0,  4.0,  5.0, 6.0  ]
        dict[ "SLDK3_T1_100K_S2_L007" ] = [ 7.0,  8.0,  9.0, 10.0, 11.0, 12.0 ]
        dict[ "SLDK3_T1_10K_S3_L007"  ] = [ 13.0, 14.0, 15.0, 16.0, 17.0, 18.0 ]
        dict[ "SLDK3_T1_1K_S4_L007"   ] = [ 19.0, 20.0, 21.0, 22.0, 23.0, 24.0 ]

        retDF = pd.DataFrame(dict)

        # we expect the gene_id to be the data frame index
        retDF.set_index("gene_id", inplace=True)

        self.logger.info(f'Count Data:\n{retDF}')
        self.logger.info(f'END')
        return retDF

    ################################################################################
    def _getMetaDataDF(self) -> pd.DataFrame :
        '''
        Get tempus, "complete seq" test meta data used to create the signature matrix
        '''
        dict = {}
        dict[ "sample_id" ] = ["SLDK3_T1_100_S1_L007", "SLDK3_T1_100K_S2_L007", "SLDK3_T1_10K_S3_L007", "SLDK3_T1_1K_S4_L007"]
        dict[ "category"  ] = ["Control", "UD", "Control", "UD"]
        dict[ "gender"    ] = ["M", "F", "M", "F"]

        retDF = pd.DataFrame(dict)

        # do not set sample_id as the index
        # retDF.set_index("sample_id", inplace=True)

        self.logger.info(f'Meta Data:\n{retDF}')
        return retDF
 
    ################################################################################
    def testCreateSignatureMatrix(self):
        '''
        TODO
        '''
        logger = logging.getLogger(__name__)
        logger.info(f'BEGIN')

        countData = self._getNormalizedCountData()

        # drop the gene_name and biotype columns
        countData = countData.drop(columns=["gene_name", "biotype"])

        metaDataDF = self._getMetaDataDF()

        retDF = createSignatureMatrix(
            genesOfInterest = ["X7D_LINE", "Zaphod"],
            normalizedCountDF = countData,
            metaDF = metaDataDF,
            useMedian = False
        )

        logger.info(f'retDF:\n{retDF}')

        expectedDF = pd.DataFrame(
            {'Control': {'X7D_LINE': 7.0, 'Zaphod': 10.0}, 
            'UD': {'X7D_LINE': 13.0, 'Zaphod': 16.0}
            } )

        # cibersort expects the index name to be 'name'
        # expectedDF.index.name = 'gene_id'
        expectedDF.index.name = 'name'

        pd.testing.assert_frame_equal(expectedDF, retDF)

        logger.info(f'END')

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
