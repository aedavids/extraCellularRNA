#
# testSignaStureGenes.py
# unit test for SignatureGenes.py 
# and pipeline.dataFactory.utilities.runSelectGenesOfInterest
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
#

import os

# from deconvolutionAnalysis.python.pipeline.dataFactory.test.signatureGeneConfigurationTest import SignatureGeneConfigurationTest
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
from pipeline.dataFactory.test.signatureGeneConfigurationTest import SignatureGeneConfigurationTest
from pipeline.dataFactory.test.testUtilities import _get1vsAllResults

from pipeline.dataFactory.driver import runSelectGenesOfInterest
import pprint as pp
import shutil
import unittest

################################################################################
class TestSignatureGenes(unittest.TestCase):
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
  
    # localCacheDir= "./data/tmp"
    localCacheDir= relativeRootPath.joinpath("data/tmp")

    expectedCandidateGenesDir = relativeRootPath.joinpath("data/testSignatureGenes/expectedCandidateGenes/testDataSet-design-tilda-gender-category-padj-001-lfc-20-n-3")
    expectedCandidateGenes = [
        str(expectedCandidateGenesDir.joinpath("UVM_vs_all.results")),
        str(expectedCandidateGenesDir.joinpath("Vagina_vs_all.results")),
        str(expectedCandidateGenesDir.joinpath("Whole_Blood_vs_all.results"))
    ]
   
    ################################################################################
    def setUp(self):
        '''
        runs before every test. Use to put system into a know state before each test
        '''
        self.logger.info("BEGIN")

        if os.path.exists(self.localCacheDir):
            self.logger.info(f"removing old localCacheDir  : {self.localCacheDir}")
            shutil.rmtree(self.localCacheDir)
        
        # outputFile = "data/geneSignatureProfiles/best/ciberSortInput/signatureGenes.tsv"
        # if os.path.exists(outputFile):
        #     self.logger.info("removing old output file : {outputFile}")
        #     os.remove(outputFile)       

        self.logger.info("END\n")

    ################################################################################
    @unittest.skip("skip template. it is not a real test")
    def testTemplate(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n")            

#    ################################################################################   
#     def _get1vsAllResults(self) -> list[str]:
#             '''
#             these are example of DESEq results file produced by 1vsAll.
#             they have multiline headers and thousands of rows (they have not been filtered)
#             '''
#             retList = [
#                 str(self.relativeRootPath.joinpath("data/testSignatureGenes/1vsAll/UVM_vs_all.results")),
#                 str(self.relativeRootPath.joinpath("data/testSignatureGenes/1vsAll/Vagina_vs_all.results")),
#                 str(self.relativeRootPath.joinpath("data/testSignatureGenes/1vsAll/Whole_Blood_vs_all.results")),
#             ]
#             return retList
            
   ################################################################################
    def testDesignPattern(self):
        '''
        test the basic design pattern

        this unit test demonstrates how a typical python driver script would be implemented to
        selected genes for a specific hypothesis
        '''
        self.logger.info("BEGIN")

        tmsg = "this is a test to see if circular imports cause problem in python"
        sgc = SignatureGeneConfigurationTest(msg=tmsg, localCacheDir=str(self.localCacheDir) )
        # sgc = SignatureGeneConfiguration(padjThreshold=0.01)

        candidateSignatureFileList = _get1vsAllResults(self.relativeRootPath)
        (selectedDict, outFileList) = runSelectGenesOfInterest(
                                        signatureGeneConfig=sgc, 
                                        candidateSignatureFileList=candidateSignatureFileList
                                        )
        tmpDir = str(self.relativeRootPath.joinpath("data/tmp/testDataSet-design-tilda-gender-category-padj-001-lfc-20-n-3"))
        expectedFileList = [
            pl.Path(tmpDir).joinpath('UVM_vs_all.results'),
            pl.Path(tmpDir).joinpath('Vagina_vs_all.results'),
            pl.Path(tmpDir).joinpath('Whole_Blood_vs_all.results')
            ]
        
        self.assertEqual( expectedFileList, outFileList)

        self.logger.info(f'selectedDict.keys() : {selectedDict.keys()}')
        expectedKeys = ['UVM_vs_all.results', 'Vagina_vs_all.results', 'Whole_Blood_vs_all.results']
        self.assertEqual(expectedKeys, list(selectedDict.keys()))
        
        expectedUVMDF = pd.read_csv(self.expectedCandidateGenes[0])
        giUVMDF = pd.read_csv(tmpDir + '/UVM_vs_all.results')
        pd.testing.assert_frame_equal(expectedUVMDF, giUVMDF)

        expectedVaginaDF = pd.read_csv(self.expectedCandidateGenes[1])
        giVDF = pd.read_csv(tmpDir + '/Vagina_vs_all.results')
        pd.testing.assert_frame_equal(expectedVaginaDF, giVDF)

        expectedWhole_Blood = pd.read_csv(self.expectedCandidateGenes[2])
        giWBDF = pd.read_csv(tmpDir + '/Whole_Blood_vs_all.results')
        pd.testing.assert_frame_equal(expectedWhole_Blood,giWBDF )

        self.logger.info(f'selectedDict : {pp.pformat(selectedDict)}')
        self.logger.info(f'outFileList : {pp.pformat(outFileList)}')

        self.logger.info("BEGIN")
      

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
