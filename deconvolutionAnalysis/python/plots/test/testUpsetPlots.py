
# testUpsetPlots.py
# unit test for upsetPlots.py 
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

import os
# AEDWIP_ENV = os.environ['AEDWIP_ENV']
# print(f'AEDWIP AEDWIP_ENV : {AEDWIP_ENV}')

ORIG_PYTHONPATH = os.environ['PYTHONPATH']
print(f'PYTHONPATH : {ORIG_PYTHONPATH}\n')
print(f'FILE: {__file__}')
print(f'PWD: {os.getcwd()}')

# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging

from   matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pathlib as pl
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration
from pipeline.dataFactory.driver import runSelectGenesOfInterest
from pipeline.dataFactory.test.testSignatureGenes import SignatureGeneConfigurationTest  
from pipeline.dataFactory.upsetPlotDataFactory import UpsetPlotDataFactory
from plots.upsetPlots import UpsetPlot
import pprint as pp
import shutil
import unittest
# import upsetplot as upsp

################################################################################
class TestUpsetPlots(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger(__name__)

    # python -m unittest discover .
    # we need to be able to find test file relative to module location
    # os.getcwd() is often the top of our source code tree
    # https://stackoverflow.com/a/57836145/4586180
    relativeRootPath = pl.Path(os.path.dirname(os.path.realpath(__file__)))

    # localCacheDir= "./data/tmp"
    localCacheDir= str(relativeRootPath.joinpath("data/tmp"))

    ################################################################################
    def setUp(self):
        '''
        runs before every test. Use to put system into a know state before each test
        '''
        self.logger.info("BEGIN")

        if os.path.exists(self.localCacheDir):
            self.logger.info(f"removing old localCacheDir  : {self.localCacheDir}")
            shutil.rmtree(self.localCacheDir)
        
        self.logger.info("END\n")

    ################################################################################
    @unittest.skip("skip template. it is not a real test")
    def testTemplate(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n")            

    ################################################################################
    def testPlot(self):
        self.logger.info("BEGIN")

        tmsg = 'plot test msg'
        title="plotTest Title"
        sgc, selectedGeneSetsDict = self._getConfiguration( tmsg, title)    
     
        upsetPlot = self._getUpsetPlot(sgc, selectedGeneSetsDict) 

        fig, pltDict = upsetPlot.configurablePlot(show_counts=True)
        plotPath = upsetPlot.savePlot( outdir=self.localCacheDir,
                                      extraFileNameParts = 'p1=2,max=3')
        self.logger.info(f'plotPath : {plotPath}')

        # test plot is correct
        expectedPngPathStr = self.relativeRootPath.joinpath('expectedUpsetPlot.png')
        self.logger.info(f"expectedPngPathStr: {expectedPngPathStr}")
        cmd = f'diff {expectedPngPathStr} {plotPath}'
        self.logger.info(f"cmd : {cmd}")
        exitStatus = os.system(cmd)
        self.logger.info(f'exitStatus : {exitStatus}')
        self.assertEqual(exitStatus, 0, 'ERROR: plot has chanaged')


        intersectionDict = upsetPlot.findIntersectionElements()
        self.logger.info(f'intersectionDict : \n{pp.pformat(intersectionDict)}')
 
        self.logger.info("END\n") 

    ################################################################################
    def testFindDegrees(self):
        '''
        TODO
        '''
        self.logger.info("BEGIN")

        tmsg = 'plot test msg'
        title="plotTest Title"
        sgc, selectedGeneSetsDict = self._getConfiguration( tmsg, title)    

        upsetPlot = self._getUpsetPlot(sgc, selectedGeneSetsDict) 

        intersectionDict = upsetPlot.findIntersectionElements()
        self.logger.info(f'intersectionDict : \n{pp.pformat(intersectionDict, indent=4, sort_dicts=True)}')

        upsetPlot.saveInteresection(self.localCacheDir, intersectionDict)

        degrees = upsetPlot.findDegrees(intersectionDict)
        self.logger.info(f'degrees: {degrees}')
        self.assertListEqual([1,2], degrees)

        for degree in degrees:
            min_degree = degree
            max_degree = degree
            fig, pltDict = upsetPlot.configurablePlot(show_counts=True, min_degree=min_degree, max_degree=max_degree)
            vargs = f"min_degree={min_degree},_max_degree={max_degree}"
            plotPath = upsetPlot.savePlot( outdir=self.localCacheDir,
                                    extraFileNameParts=vargs)
            self.logger.info(f'plotPath : {plotPath}')
            
        self.logger.info("END")

    ################################################################################
    def testIntersection(self):
        '''
        The upsetPlot shows set level information of intersection. We want to know what the
        actual elements in the intersection are
        '''
        self.logger.info("BEGIN")

        tmsg = 'plot test msg'
        title = "intersection" 
        sgc, selectedGeneSetsDict = self._getConfiguration(tmsg, title)

        upsetPlot = self._getUpsetPlot(sgc, selectedGeneSetsDict) 

        intersectionDict = upsetPlot.findIntersectionElements()   

        print("\n\n+++++++++++++")
        print(f'intersectionDict:\n{pp.pformat(intersectionDict)}')
        
        expectedIntersectionDict = self._getExpectectedIntersectionDict()
        self.assertDictEqual(expectedIntersectionDict, intersectionDict, msg="ERROR intersectionDict does not match expected")

        self.logger.info("END")

    ###############################################################################  
    def _getConfiguration(self, msg, title) -> (SignatureGeneConfigurationTest, dict):
        sgc = SignatureGeneConfigurationTest(msg=msg, 
                                            localCacheDir=self.localCacheDir, 
                                            title=title)

        candidateSignatureFileList = self._get1vsAllResults()

        selectedGeneSetsDict, outFileList = runSelectGenesOfInterest(
                                        signatureGeneConfig=sgc, 
                                        candidateSignatureFileList=candidateSignatureFileList 
                                        )
        self.logger.info(f'selectedGeneSetsDict:\n{pp.pformat(selectedGeneSetsDict)}')

        return (sgc, selectedGeneSetsDict)

    ###############################################################################  
    def _getExpectectedIntersectionDict(self) :
            # expectedIntersectionDict = {
            #                     'UVM': ['UVM_AAA', 'UVM_V_W.1'],
            #                     'UVM_XXX_Vagina': ['UVM_V'],
            #                     'Vagina': ['V_BBB'],
            #                     'Vagina_XXX_Whole_Blood': ['UVM_V_W'],
            #                     'Whole_Blood': ['W_CCC', 'V_W']
            # }

            expectedIntersectionDict = {
                ('Whole_Blood',): ['W_CCC', 'V_W'],
                ('Vagina',): ['V_BBB'],
                ('Vagina', 'Whole_Blood'): ['UVM_V_W'],
                ('UVM',): ['UVM_AAA', 'UVM_V_W.1'],
                ('UVM', 'Vagina'): ['UVM_V']                                
            }

            return expectedIntersectionDict

    ###############################################################################   
    def _get1vsAllResults(self) -> list[str]:
            '''
            these are example of DESEq results file produced by 1vsAll.
            they have multiline headers and thousands of rows (they have not been filtered)
            '''            
            retList = [
                str( self.relativeRootPath.joinpath("../../pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/UVM_vs_all.results") ),
                str( self.relativeRootPath.joinpath("../../pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/Vagina_vs_all.results") ),
                str( self.relativeRootPath.joinpath("../../pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/Whole_Blood_vs_all.results") ),
            ]
            return retList

    ###############################################################################   
    def _getUpsetPlot(self, 
                      sgc : SignatureGeneConfigurationTest,
                      selectedGeneSetsDict : dict):
        upsetFactory = UpsetPlotDataFactory()
        upsetPlotDataDF, geneSetsDict = upsetFactory.createUpsetPlotDataFromSetDict(selectedGeneSetsDict)
        print("\n\n\n###################")

        upsetPlot = UpsetPlot(signatureGeneConfig=sgc, 
                    geneSetsDict=geneSetsDict, 
                    geneSetsUpsetPlotData=upsetPlotDataDF)   
        
        return upsetPlot
    
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
