'''
Created on Jun 19, 2020

@author: andrewdavidson
'''
import logging
from   setupLogging import setupLogging
import unittest


################################################################################
class TestLoggingConfi( unittest.TestCase ):
    configFilePath = setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)
    logger.info("using logging configuration file:{}".format(configFilePath)) 
    
    ################################################################################
    def setUp( self ):
        pass

    ################################################################################
    def tearDown( self ):
        pass

    ################################################################################
    def testTemplate( self ):
        self.logger.info("BEGIN")
        
        self.logger.debug("AEDWIP debug")
        self.logger.info("AEDWIP info")
        self.logger.warning("AEDWIP warn:{}".format(42))
        self.logger.error("AEDWIP error:{}".format(45))
        self.logger.info("END\n")

    ################################################################################
    def testConfig( self ):
        self.logger.info("BEGIN")
        
        self.logger.info("END\n")


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
