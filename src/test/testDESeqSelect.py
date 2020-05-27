'''
Created on May 26, 2020

@author: andrewdavidson
'''

from kimLabDEQ.DESeqSelect import DESeqSelect
import logging
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
            
        self.logger.info("END")            
            
    ################################################################################
    def testName(self):
        self.logger.info("BEGIN")
        
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.logger.info('created temporary directory:{}'.format(tmpdirname))
            adjPValue = 0.7
            logFoldChange = 0.03
            s = DESeqSelect("data/testGencode.de.deq.csv", tmpdirname, adjPValue, logFoldChange)
            #s = DESeqSelect("data/testGencode.de.deq.csv", 'aedwip', adjPValue, logFoldChange)
            
            outputFileName = s.getOutputFileName()
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
        
            
################################################################################
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()