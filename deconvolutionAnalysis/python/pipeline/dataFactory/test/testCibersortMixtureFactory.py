#
# testCibersortSignatureMatrixFactory.py
# unit test for TestCibersortSignatureMatrixFactory.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/createCibersortMixtureMatrix.ipynb
#

import os
import pandas as pd
import pathlib as pl
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
from pipeline.dataFactory.cibersortMixtureMatrixFactory import CibersortMixtureFactory
import shutil
import unittest

################################################################################
class TestCibersortMixtureMatrixFactory(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger(os.path.basename(__file__))

    # https://stackoverflow.com/a/14493895/4586180
    #maxDiff = None

    # python -m unittest discover .
    # we need to be able to find test file relative to module location
    # os.getcwd() is often the top of our source code tree
    # https://stackoverflow.com/a/57836145/4586180
    relativeRootPath = pl.Path(os.path.dirname(os.path.realpath(__file__)))
 
    localCacheDir= relativeRootPath.joinpath("data/tmp")
    outdir = relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput")
    outputFiles = [
                    relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/testCibersortMixtureFactory_mixture.txt"),
                    relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/testCibersortMixtureFactory_expectedFractions.txt"),
                    relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/testCibersortMixtureFactory_randomizedMixture.txt")
        ]
    prefix = "testCibersortMixtureFactory"

    ################################################################################
    def setUp(self):
        '''
        runs before every test. Use to put system into a know state before each test
        '''
        self.logger.info("BEGIN")

        if os.path.exists(self.localCacheDir):
            self.logger.info(f"removing old localCacheDir  : {self.localCacheDir}")
            shutil.rmtree(self.localCacheDir)       

        for outputFile in self.outputFiles:
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

        # use the testCibersortMatrixFactory expected output
        signatueGeneFilePath = self.relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/expectedSignatureGenes.tsv")
        groupByGeneCountFilePath =  self.relativeRootPath.joinpath("data/trainGroupby.csv")
        colDataFilePath = self.relativeRootPath.joinpath("data/trainColData.csv")
        scalingFactorsPath = self.relativeRootPath.joinpath("data/1vsAll/estimatedSizeFactors.csv")
        csmf = CibersortMixtureFactory(
                    str(signatueGeneFilePath),
                    str(groupByGeneCountFilePath), 
                    str(colDataFilePath),
                    str(scalingFactorsPath),
                    str(self.localCacheDir)
                )
        
        mixturePath, expectedFractionsPath = csmf.saveMixtureAndExpectedFractions(outDir=self.outdir, prefixStr=self.prefix)
        self.assertEqual(self.outputFiles[0], mixturePath)
        self.assertEqual(self.outputFiles[1], expectedFractionsPath)

        expectedPath = self.relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/expectedTestCibersortMixtureFactory_expectedFractions.txt")
        expectedFractionDF = pd.read_csv(expectedPath, sep='\t')
        pd.testing.assert_frame_equal(expectedFractionDF, csmf.getExpectedFractionsDF())

        # aedwip add assert for mixture

        expectedMixturePath =  self.relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/expectedTestCibersortMixtureFactory_mixture.txt")
        expectedMixtureDF = pd.read_csv(expectedMixturePath, sep='\t')
        labeledMixtureDF = csmf.getLabeledMixtureDF()
        self.logger.info(f'labeledMixtureDF\n{labeledMixtureDF}')

        mixtureDF = csmf._convertToMixtureToCiberSortFmt( csmf.labeledMixtureDF )
        # to_csv does not save column name attribute
        mixtureDF.columns.name = None
        self.logger.info(f'mixtureDF\n{mixtureDF}')
        pd.testing.assert_frame_equal(expectedMixtureDF, mixtureDF)

        # this is not a great test
        # the values in each row are the same. so randomize does not change anythin
        randomizedDF = csmf.randomizeMixture(seed=42)
        randomizedPath = csmf.saveRandomizedMixture(outDir="data/geneSignatureProfiles/best/ciberSortInput", 
                                                    randomizedDF=randomizedDF,
                                                    prefixStr=self.prefix)
        
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
