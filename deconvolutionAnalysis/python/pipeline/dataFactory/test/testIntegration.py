#
# testIntegration.py
# end to end test we can use to debug slurm scripts
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

# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging

import numpy as np
import pandas as pd
import pathlib as pl

from pipeline.dataFactory.cibersortMixtureMatrixFactory import CibersortMixtureFactory
from pipeline.dataFactory.cibersortSignatureMatrixFactory import CibersortSignatureMatrixFactory
from pipeline.dataFactory.driver import runSelectGenesOfInterest
from pipeline.dataFactory.test.signatureGeneConfigurationTest import SignatureGeneConfigurationTest
from pipeline.dataFactory.test.testUtilities import _get1vsAllResults

from pipeline.dataFactory.upsetPlotDataFactory import UpsetPlotDataFactory
from plots.upsetPlots import UpsetPlot

import pprint as pp
import shutil
import unittest

################################################################################
class TestIntegration(unittest.TestCase):
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
  
    # localCacheDir= "./data/tmp"
    localCacheDir= relativeRootPath.joinpath("data/testIntegration/tmp")

    # The files produced by runSelectGenesOfInterest() should be found here
    filteredDESeqResultsDir = localCacheDir.joinpath('testDataSet-design-tilda-gender-category-padj-001-lfc-20-n-3')
   
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
    def testPart1(self):
        '''
        prepare data for CIBERSORTx
            1. select genes of interest
            2. create upset plots and find interesection elements
            3. create CIBERSORTx signatureGenes matrix
            4. create CIBERSORTx mixture matrix
        '''
        self.logger.info("BEGIN")

        #
        # step 1
        #
        self.logger.info("step 1 select genes of interest")
        tmsg = "integration unit test"
        sgc = SignatureGeneConfigurationTest(msg=tmsg, localCacheDir=str(self.localCacheDir) )
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

        self.logger.info(f'outFileList:\n{pp.pformat(outFileList)}')    
        self.logger.info(f'selectedGeneSetsDict:\n{pp.pformat(selectedGeneSetsDict)}')    

        #
        # step 2
        #
        self.logger.info('step 2 create plots and find intersection elements')
        upsetFactory = UpsetPlotDataFactory()
        upsetPlotDataDF, geneSetsDict = upsetFactory.createUpsetPlotDataFromSetDict(selectedGeneSetsDict)
        self.logger.info(f'upsetPlotDataDF:\n{pp.pformat(upsetPlotDataDF)}')

        upsetPlot = UpsetPlot(signatureGeneConfig=sgc, 
                        geneSetsDict=geneSetsDict, 
                        geneSetsUpsetPlotData=upsetPlotDataDF)
        
     

        fig, pltDict = upsetPlot.configurablePlot(show_counts=True)
        upsetPlotOutDir = os.path.join(outDir, "upsetPlot.out")
        plotPath = upsetPlot.savePlot( outdir=upsetPlotOutDir,
                                      extraFileNameParts = 'p1=2,max=3')
        self.logger.info(f'plotPath : {plotPath}')

        intersectionDict = upsetPlot.findIntersectionElements()  

        intersectionDictPath = upsetPlot.saveInteresection(upsetPlotOutDir, intersectionDict) 
        self.logger.info(f'intersectionDictPath : {intersectionDictPath}')

        expectedIntersectionDict = self._expectedIntersectionDict()
        self.assertDictEqual(expectedIntersectionDict, intersectionDict, msg="ERROR intersectionDict does not match expected")


        #
        # step 3
        # create data into CIBERSORTx signatureGenes matrix
        #
        self.logger.info('step 3: create CIBERSORTx signatureGenes matrix')

        self.logger.info(f'sgc.getLocalCachedDir():      {sgc.getLocalCachedDir()}')
        self.logger.info(f'self.filteredDESeqResultsDir: {self.filteredDESeqResultsDir}')
        emsg = "ERROR the filtered DESeq results file are not in the expected directory"
        self.assertEqual(str(self.filteredDESeqResultsDir), sgc.getLocalCachedDir(), emsg )

        countPath = self.relativeRootPath.joinpath("data/testIntegration/geneCounts.csv")
        colDataPath = self.relativeRootPath.joinpath("data/testIntegration/colData.csv")
        esfPath = self.relativeRootPath.joinpath("data/testSignatureGenes/1vsAll/estimatedSizeFactors.csv")
        csmf = CibersortSignatureMatrixFactory(
            geneSignatureProfilesDataRootDir = outDir, #str( self.filteredDESeqResultsDir ),
            #oneVsAllDataDir = "data/1vsAll", 
            groupByGeneCountFilePath = str(countPath),
            colDataFilePath = str(colDataPath), 
            estimatedScalingFactorsFilePath = str(esfPath), 
            localCacheDir= outDir, #str(self.localCacheDir),
            outdir = "ciberSortInput", 
            testSize = None, 
            verbose = False
        )
            
        signatureGeneDF = csmf.getCiberSortSignatueDF()

        self.logger.info("signatureGeneDF.shape:{}".format(signatureGeneDF.shape)) 
        self.logger.info( f'signatureGeneDF\n {signatureGeneDF}' )
        self.assertEqual( (7,3), signatureGeneDF.shape)
        pathToSignatureFile = csmf.save()

        pathToSignatureFile = self.relativeRootPath.joinpath(pathToSignatureFile)
        self.logger.info(f'pathToSignatureFile : {pathToSignatureFile}')

        expectedSignatureFilePath = self.localCacheDir.joinpath("testDataSet-design-tilda-gender-category-padj-001-lfc-20-n-3/ciberSortInput/signatureGenes.tsv")
        self.assertEqual(expectedSignatureFilePath, pathToSignatureFile) 

        expectedSignatureDF = self._getExpectedSignatureDF()
        self.logger.info(f'expectedSignature Data Frame\n{expectedSignatureDF}')
        self.logger.info(f'signatureGeneDF Data Frame\n{signatureGeneDF}')
        pd.testing.assert_frame_equal(expectedSignatureDF, signatureGeneDF)        

        #
        # step 4
        # create CIBERSORTx mixture matrix
        #
        self.logger.info('step 4: create CIBERSORTx mixture matrix')
    
        csmf = CibersortMixtureFactory(
                    str(pathToSignatureFile),
                    str(countPath), #str(groupByGeneCountFilePath), 
                    str(colDataPath), #str(colDataFilePath),
                    str(esfPath), #str(scalingFactorsPath),
                    str(self.localCacheDir)
                )
        
        prefix = "testCibersortMixtureFactory"
        saveDir = os.path.join(outDir, "ciberSortInput")
        mixturePath, expectedFractionsPath = csmf.saveMixtureAndExpectedFractions(outDir=saveDir, prefixStr=prefix)
        self.logger.info(f'mixturePath : {mixturePath}')
        self.logger.info(f'expectedFractionsPath : {expectedFractionsPath}')

        randomizedMixtureDF = csmf.randomizeMixture(seed=42)
        randomizedPath = csmf.saveRandomizedMixture(outDir=saveDir, randomizedDF=randomizedMixtureDF, prefixStr=prefix)
        self.logger.info(f'randomizedPath : {randomizedPath}')

        self.logger.info("END\n")            

################################################################################
    def _getExpectedSignatureDF(self):
        '''
        knock outs are correct

        looks like counts are in correct ball park. TODO verify
        '''
      
        retDF = pd.DataFrame(
                    {
                        'UVM': {'UVM_AAA': 1.6,
                                'UVM_V': 1.6,
                                'UVM_V_W': 0.0,
                                'UVM_V_W.1': 1.6,
                                'V_BBB': 0.0,
                                'V_W': 0.0,
                                'W_CCC': 0.0},

                    'Vagina': {'UVM_AAA': 0.0,
                                'UVM_V': 11.1,
                                'UVM_V_W': 11.1,
                                'UVM_V_W.1': 0.0,
                                'V_BBB': 11.1,
                                'V_W': 0.0,
                                'W_CCC': 0.0},

                    'Whole_Blood': {'UVM_AAA': 0.0,
                                    'UVM_V': 0.0,
                                    'UVM_V_W': 29.0,
                                    'UVM_V_W.1': 0.0,
                                    'V_BBB': 0.0,
                                    'V_W': 29.0,
                                    'W_CCC': 29.0}
                        }
                )

        retDF.index.name = "name"
        return retDF
    
################################################################################
    def _expectedIntersectionDict(self):
         expectedIntersectionDict = {
                                'UVM': ['UVM_AAA', 'UVM_V_W.1'],
                                'UVM_XXX_Vagina': ['UVM_V'],
                                'Vagina': ['V_BBB'],
                                'Vagina_XXX_Whole_Blood': ['UVM_V_W'],
                                'Whole_Blood': ['W_CCC', 'V_W']
            }
         
         return expectedIntersectionDict
    
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
