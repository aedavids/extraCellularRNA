#
# testRankSignatureGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

import os
# AEDWIP_ENV = os.environ['AEDWIP_ENV']
# print(f'AEDWIP AEDWIP_ENV : {AEDWIP_ENV}')

ORIG_PYTHONPATH = os.environ['PYTHONPATH']
print(f'PYTHONPATH : {ORIG_PYTHONPATH}')
print(f'FILE: {__file__}')
print(f'PWD: {os.getcwd()}')

import ast

# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging

# import numpy as np
import pandas as pd
import pathlib as pl
import pickle
import pprint as pp
import shutil
import unittest

from analysis.rankSignatureGeneConfig import RankSignatureGeneConfig
from pipeline.dataFactory.driver import runSelectGenesOfInterest


################################################################################
class TestByRankSignatureGeneConfig(unittest.TestCase):
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
    localCacheDir= relativeRootPath.joinpath("data/testRankSignatureGeneConfig/tmp")

    # The files produced by runSelectGenesOfInterest() should be found here
    filteredDESeqResultsDir = localCacheDir.joinpath('testDataSet-design-tilda-gender-category-padj-001-lfc-20-n-3')

    # expectedCandidateGenesDir = relativeRootPath.joinpath("data/testEnrichSignatureGeneConfig/1vsAllExpected")
    # expectedCandidateGenes = [
    #     str(expectedCandidateGenesDir.joinpath("UVM_vs_all.results")),
    #     str(expectedCandidateGenesDir.joinpath("Vagina_vs_all.results")),
    #     str(expectedCandidateGenesDir.joinpath("Whole_Blood_vs_all.results"))
    # ]
   
    ################################################################################
    def setUp(self):
        '''
        runs before every test. Use to put system into a know state before each test
        '''
        self.logger.info("BEGIN")

        if os.path.exists(self.localCacheDir):
            self.logger.info(f"removing old localCacheDir  : {self.localCacheDir}")
            shutil.rmtree(self.localCacheDir)
        
        signatureMatrixOutputFile = self.relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/signatureGenes.tsv")
        if os.path.exists(signatureMatrixOutputFile):
            self.logger.info(f"removing old output file : {signatureMatrixOutputFile}")
            os.remove(signatureMatrixOutputFile)           

        self.logger.info("END\n")

    ################################################################################
    @unittest.skip("skip template. it is not a real test")
    def testTemplate(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n")            

  ################################################################################
    def testRankSignatureGeneConfig(self):
        '''
        deconvolutionAnalysis/python/analysis/test/data/testEnrichSignatureGeneConfig/intersection.dict
        
            ['UVM_AAA', 'UVM_V, 'UVM_V_W', 'UVM_V_W.1', 'V_BBB', 'V_W',  'W_CCC', ],
        wb        YYY                 X                                X       X
        V                    X        X            YYY        X        YYY     YYY
        UVM        X         X                     X

        intersections
        Whole_Blood and UVM do not have any unique genes
        {
            ('Vagina',)               : ['V_BBB'],
            ('Vagina', 'Whole_Blood') : ['UVM_V_W', 'V_W', 'W_CCC'],
            ('UVM', 'Vagina')         : ['UVM_V', 'UVM_V_W.1'],
            ('UVM', 'Whole_Blood')    : ['UVM_AAA']
        }      
        '''
        
        self.logger.info("BEGIN")
        tmsg = "enrich unit test"
        # sgc = SignatureGeneConfigurationTest(msg=tmsg, localCacheDir=str(self.localCacheDir) )
        dataSetName = "testDataSet"
        design = "~gender + category"
        padjThreshold =0.01
        lfcThreshold = 2.0
        n = 8
        localCacheRootPath =str(self.localCacheDir)
        title ="test enrichment"
        deseqResultsDir = str( self.relativeRootPath.joinpath("data/testEnrichSignatureGeneConfig/1vsAll") )
       
        sgc = RankSignatureGeneConfig(
                dataSetName, 
                design, 
                padjThreshold, 
                lfcThreshold, 
                n, 
                localCacheRootPath, 
                title,
                deseqResultsDir,
        )

        outDir = sgc.getLocalCachedDir()
        self.logger.info(f'outDir : {outDir}')

        # runSelectGenesOfInterest
        # filters/selects the genes of interest from each of the 1vsAll DESeq  results file
        # the selectedGeneSetsDict contains a key/value for each file
        # the outFileList is a file with a list of path to filtered DESeq results file.
        # ie only DESeq results we are interested
        candidateSignatureFileList = _get1vsAllResults(self.relativeRootPath)
        (selectedGeneSetsDict, outFileList) = runSelectGenesOfInterest(
                                        signatureGeneConfig=sgc, 
                                        candidateSignatureFileList=candidateSignatureFileList 
                                        ) 
        
        self.logger.info(f'outFileList :\n{outFileList}')
        self.logger.info(f'selectedGeneSetsDict :\n{pp.pformat(selectedGeneSetsDict)}')

        expectedRankedDict = {
            'Whole_Blood': ['V_W', 'W_CCC', 'UVM_V_W', 'ALAS2', 'CXCR1', 'HBG2', 'HBD', 'AC104389.6'], 
            'Vagina': ['C1orf61', 'ARPP21', 'OLIG1', 'OLIG2', 'SILC1', 'UGT2A3', 'UVM_V', 'V_BBB'], 
            'UVM': ['AC114812.1', 'UVM_V_W.1', 'UVM_AAA', 'LINC02612', 'AC090150.1', 'ANO1-AS1', 'PCSEAT']
            }

        expectedRankedDict = {
            'Whole_Blood': ['V_W', 'W_CCC', 'UVM_V_W', 'ALAS2', 'CXCR1', 'HBG2', 'HBD', 'AC104389.6'], 
            'Vagina': ['C1orf61', 'ARPP21', 'OLIG1', 'OLIG2', 'SILC1', 'UGT2A3', 'UVM_V', 'V_BBB'], 
            'UVM': ['AC114812.1', 'UVM_V_W.1', 'UVM_AAA', 'LINC02612', 'AC090150.1', 'ANO1-AS1', 'PCSEAT', 'AC244102.4']
        }
        self.logger.info(f'sgc.rankedDict : \n {sgc.rankedDict}')
        self.assertDictEqual(expectedRankedDict, sgc.rankedDict)

        expectedRankedPaths = _getExpected1vsAllPaths(self.relativeRootPath)
        eUVMPath = expectedRankedPaths[0]
        eVaginaPath = expectedRankedPaths[1]
        eWholeBloodPath = expectedRankedPaths[2]

        uvmDF = selectedGeneSetsDict['UVM_vs_all.results']
        vaginaDF = selectedGeneSetsDict['Vagina_vs_all.results']
        wholeBloodDF = selectedGeneSetsDict['Whole_Blood_vs_all.results']

        # uvmDF.to_csv(eUVMPath, index=False)
        # vaginaDF.to_csv(eVaginaPath, index=False)
        # wholeBloodDF.to_csv(eWholeBloodPath, index=False)
        
        expectedUVMDF = pd.read_csv(eUVMPath)
        uvmDF.reset_index(drop=True, inplace=True)
        pd.testing.assert_frame_equal(expectedUVMDF, uvmDF)

        expectedVaginaDF = pd.read_csv(eVaginaPath)
        vaginaDF.reset_index(drop=True, inplace=True)
        pd.testing.assert_frame_equal(expectedVaginaDF, vaginaDF)

        expectedWholeBloodDF = pd.read_csv(eWholeBloodPath)
        wholeBloodDF.reset_index(drop=True, inplace=True)
        pd.testing.assert_frame_equal(expectedWholeBloodDF, wholeBloodDF)


        self.logger.info("END")   

################################################################################   
def _get1vsAllResults(relativeRootPath : pl.Path) -> list[str]:
    '''
    these are example of DESEq results file produced by 1vsAll.
    they have multiline headers and thousands of rows (they have not been filtered)
    '''
    retList = [
        str(relativeRootPath.joinpath("data/testEnrichSignatureGeneConfig/1vsAll/UVM_vs_all.results")),
        str(relativeRootPath.joinpath("data/testEnrichSignatureGeneConfig/1vsAll/Vagina_vs_all.results")),
        str(relativeRootPath.joinpath("data/testEnrichSignatureGeneConfig/1vsAll/Whole_Blood_vs_all.results")),
    ]
    return retList

################################################################################   
def _getExpected1vsAllPaths(relativeRootPath : pl.Path) -> list[str]:
        '''
        these are example of DESEq results file produced by 1vsAll.
        they have multiline headers and thousands of rows (they have not been filtered)
        '''
        retList = [
            str(relativeRootPath.joinpath("data/testRankSignatureGeneConfig/expected1vsAll/UVM_vs_all.results")),
            str(relativeRootPath.joinpath("data/testRankSignatureGeneConfig/expected1vsAll/Vagina_vs_all.results")),
            str(relativeRootPath.joinpath("data/testRankSignatureGeneConfig/expected1vsAll/Whole_Blood_vs_all.results")),
        ]
        return retList


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
