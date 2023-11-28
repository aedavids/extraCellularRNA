#
# testCibersortSignatureMatrixFactory.py
# unit test for TestCibersortSignatureMatrixFactory.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
#

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

#from analysis.fractions import CibersortResultsAsKWayClassifier
import numpy as np
import pandas as pd
import pathlib as pl
from pipeline.dataFactory.cibersortSignatureMatrixFactory import CibersortSignatureMatrixFactory
import shutil
import unittest

################################################################################
class TestCibersortSignatureMatrixFactory(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger(__name__)

    # https://stackoverflow.com/a/14493895/4586180
    #maxDiff = None

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
    def getExpectedSignatureDF(self) -> pd.DataFrame :
        retDF = pd.read_csv(self.expectedSignatureFilePath, sep="\t", index_col="name")
        return retDF
    
    ################################################################################
    def setUp(self):
        '''
        runs before every test. Use to put system into a know state before each test
        '''
        self.logger.info("BEGIN")

        if os.path.exists(self.localCacheDir):
            self.logger.info(f"removing old localCacheDir  : {self.localCacheDir}")
            shutil.rmtree(self.localCacheDir)
        
        outputFile = self.relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/signatureGenes.tsv")
        if os.path.exists(outputFile):
            self.logger.info(f"removing old output file : {outputFile}")
            os.remove(outputFile)       

        self.logger.info("END\n")

    ################################################################################
    @unittest.skip("skip template. it is not a real test")
    def testTemplate(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n")            
                
    ################################################################################
    def testFactory(self):
        self.logger.info("BEGIN")

        csmf = CibersortSignatureMatrixFactory(
            geneSignatureProfilesDataRootDir = str(self.relativeRootPath.joinpath("data/geneSignatureProfiles/best")),
            #oneVsAllDataDir = "data/1vsAll", 
            groupByGeneCountFilePath = str(self.relativeRootPath.joinpath("data/trainGroupby.csv")),
            colDataFilePath = str(self.relativeRootPath.joinpath("data/trainColData.csv")), 
            estimatedScalingFactorsFilePath = str(self.relativeRootPath.joinpath("data/1vsAll/estimatedSizeFactors.csv")), 
            localCacheDir= str(self.localCacheDir),
            outdir = "ciberSortInput", 
            testSize = None, 
            verbose = False
        )
            
        signatureGeneDF = csmf.getCiberSortSignatueDF()

        #
        # goofy hack to taht assert will work
        # for unknow reason 
        self.logger.info("signatureGeneDF.shape:{}".format(signatureGeneDF.shape))
        self.logger.info( f'bestSignatureGeneDF\n {signatureGeneDF}' )
        pathToSignatureFile = csmf.save()

        pathToSignatureFile = self.relativeRootPath.joinpath(pathToSignatureFile)
        self.logger.info(f'pathToSignatureFile : {pathToSignatureFile}')

        self.assertEqual(self.expectedSignatureFilePath, pathToSignatureFile)

        expectedSignatureDF = self.getExpectedSignatureDF()
        self.logger.info(f'expectedSignature Data Frame\n{expectedSignatureDF}')
        self.logger.info(f'signatureGeneDF Data Frame\n{signatureGeneDF}')
        pd.testing.assert_frame_equal(expectedSignatureDF, signatureGeneDF)

        self.logger.info("END\n") 

    ################################################################################
    '''
    make sure setup() just before every run. check log file
    '''
    def testAAASetUp(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n" )         

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
