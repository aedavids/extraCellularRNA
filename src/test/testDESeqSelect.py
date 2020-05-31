'''
Created on May 26, 2020

@author: andrewdavidson
'''

from kimLabDEQ.DESeqSelect import DESeqSelect
import logging
import numpy as np
from   setupLogging import setupLogging
import tempfile
import unittest

################################################################################
class TestDESeqSelect(unittest.TestCase):
    configFilePath = setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)
    logger.info("using logging configuration file:{}".format(configFilePath))   
    
    ################################################################################
    def testTemplate(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n")            
            
    ################################################################################
    def testSaveRows(self):
        self.logger.info("BEGIN")
        
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.logger.info('created temporary directory:{}'.format(tmpdirname))
            adjPValue = 0.7
            logFoldChange = 0.03
            s = DESeqSelect("data/testGencode.de.deq.csv")
            outputFileName = s.saveRows( tmpdirname, adjPValue, logFoldChange)
            #s = DESeqSelect("data/testGencode.de.deq.csv", 'aedwip', adjPValue, logFoldChange)
            
            self.logger.info("outputFileName:{}".format(outputFileName))
                
            self.logger.info("END")            
            
            results = []
            with open(outputFileName, "r") as resultsFH:
                for resultLine in resultsFH.readlines():
                    results.append(resultLine)
                    print(resultLine)
                
            expectedLines = [
                "\"\",\"baseMean\",\"log2FoldChange\",\"lfcSE\",\"pvalue\",\"padj\"\n",
                "\"A2ML1\",9.06198048230035,0.158921267550763,0.468595251083387,0.344503244538388,0.70745599036947\n",
                "\"AACS\",153.862914271551,0.0386190650063446,0.278544076196668,0.865165709209835,0.959407298838613\n"
                ]

            self.assertListEqual(expectedLines, results)
            
        self.logger.info("END\n")            
        
    ################################################################################
    def testReadVolcanoPlotData(self):
        self.logger.info("BEGIN")
           
        s = DESeqSelect("data/testGencode.de.deq.csv")
        geneNamesNP, xNP, yNP = s.readVolcanoPlotData()
        expectedLogFold = np.array([-0.20356742, -0.35149798, -0.10502733,  0.15892127, -0.13266999,  
                                    0.280224,   0.03861907])
        
        self.logger.info("geneNames:{}".format(geneNamesNP))
        
        self.logger.info("xNP:\n{}".format(xNP))
        self.logger.info("expectedLogFold:\n{}".format(expectedLogFold))
        
        self.logger.info("yNP:{}".format(yNP))
        
        expectedGenes = ['"7SK"' '"A1BG"' '"A2M"' '"A2ML1"' '"A4GALT"' '"AAAS"' '"AACS"']

        expectedLogAdjP = np.array([0.13261,      0.30662471,  0.12950182,  0.15030057,  0.13678097, 
                            0.21152287, 0.01799698])
        
        self.assertListEqual(expectedGenes, expectedGenes)
        self.assertTrue( np.allclose(expectedLogFold, xNP) ) 
        self.assertTrue( np.allclose(expectedLogAdjP, yNP) ) 

        self.logger.info("END\n")          
            
################################################################################
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()