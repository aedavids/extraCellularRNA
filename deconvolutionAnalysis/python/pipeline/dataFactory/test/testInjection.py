#
# testInjection.py
# Can we dynamical load a file and call a function?
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: https://stackoverflow.com/questions/301134/how-can-i-import-a-module-dynamically-given-its-name-as-string
#

import importlib
import os
# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging
import pathlib as pl
import shutil
import unittest

# from deconvolutionAnalysis.python.pipeline.dataFactory.test.signatureGeneConfigurationTest import SignatureGeneConfigurationTest
# AEDWIP_ENV = os.environ['AEDWIP_ENV']
# print(f'AEDWIP AEDWIP_ENV : {AEDWIP_ENV}')

ORIG_PYTHONPATH = os.environ['PYTHONPATH']
print(f'PYTHONPATH : {ORIG_PYTHONPATH}')
print(f'FILE: {__file__}')
print(f'PWD: {os.getcwd()}')


################################################################################
class TestInjection(unittest.TestCase):
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
    localCacheDir= relativeRootPath.joinpath("data/tmp")

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
    def testDynamicLoading(self):
        self.logger.info("BEGIN")

        try:
            name = "pipeline.dataFactory.test.signatureGeneConfigurationTest"
            myModule = importlib.import_module(name, package=None)
            sgct = myModule.SignatureGeneConfigurationTest("msg", str(self.localCacheDir))
            self.logger.info(f'sgct : {sgct}')
            retStr = sgct.testDynamicalyloadedModules()
            self.logger.info(f'sgct.testDynamicalyloadedModules() : {retStr}')
            self.assertEqual(retStr, "hello world")

        except ImportError as e:
            emsg = f'unable to dynamically import {name} exception: {e}'
            self.logger.error(emsg)
            self.assertFalse(True, emsg)


 
        self.logger.info("END\n")         



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
