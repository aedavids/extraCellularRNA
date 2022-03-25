'''
Created on Jan 19, 2022

@author: andrewdavidson
'''
# use findspark to inject pyspark into process
import findspark
import os
from bigDataDeseq.countMatrix import CountMatrix
# import bigDataDeseq

SPARK_HOME_ENV = os.getenv( "SPARK_HOME" )
if SPARK_HOME_ENV:
    # logger.info("using SPARK_HOME environment variable:{}".format(SPARK_HOME))
    findspark.init()
else:
    SPARK_HOME_HACK = "../../spark-3.1.2-bin-hadoop3.2"
    findspark.init( SPARK_HOME_HACK )

import logging
import pandas as pd
from   pyspark.sql import SparkSession
from   bigDataDeseq.setupLogging import setupLogging
import unittest

################################################################################
class TestCountMatrix(unittest.TestCase):
    spark = SparkSession\
                .builder\
                .appName( "TestEstimatedScalingFactors" )\
                .getOrCreate()
                    # .config("spark.driver.memory", "15g")     .getOrCreate()

    # this pyspark logging initialization works how ever does not produce nicely formated messages
    # TODO read how to configure output
    # https://medium.com/@lubna_22592/building-production-pyspark-jobs-5480d03fd71e
    # initialize  logger for yarn cluster logs
    # log4jLogger = sc._jvm.org.apache.log4j
    # log4jLogger = spark.sparkContext._jvm.org.apache.log4j
    # logger = log4jLogger.LogManager.getLogger(__name__)
    # logger.info("pyspark script logger initialized")    
    
    logger = logging.getLogger( __name__ )
    configFilePath = setupLogging( default_path='logging.test.ini.json' )
    logger.info( "using logging configuration file:{}".format( configFilePath ) )
    logger.info( "current working directory:{}".format( os.getcwd() ) )

    if SPARK_HOME_ENV:
        logger.info( "using SPARK_HOME environment variable:{}".format( SPARK_HOME_ENV ) )
    else:
        logger.warning( "using SPARK_HOME={}".format( SPARK_HOME_HACK ) )
        logger.warning( "set SPARK_HOME environment variable to change spark version or install location" )

    fileList = [
                    "data/ctrl_1.quant.sf",
                    "data/ctrl_2.quant.sf",
                    "data/ctrl_3.quant.sf",
                    "data/kras_1.quant.sf",
                    "data/kras_2.quant.sf",
                    "data/kras_3.quant.sf"
            ]
    
    sampleNames = [ "ctrl-1", "ctrl_2", "ctrl_3", "kras_1", "kras_2", "kras_3" ]    
    
    expectedMatrixDict = {
        'Name': {0: 'txId_1', 1: 'txId_2', 2: 'txId_3', 3: 'txId_4', 4: 'txId_5', 5: 'txId_6', 6: 'txId_7', 7: 'txId_8', 8: 'txId_9', 9: 'txId_10'},
        'ctrl-1': {0: 0.0, 1: 11.0, 2: 12.0, 3: 13.0, 4: 14.0, 5: 15.0, 6: 16.0, 7: 17.0, 8: 18.0, 9: 19.0},
        'ctrl_2': {0: 0.1, 1: 11.1, 2: 12.1, 3: 13.1, 4: 14.1, 5: 15.1, 6: 16.1, 7: 17.1, 8: 18.1, 9: 19.1},
        'ctrl_3': {0: 0.2, 1: 11.2, 2: 0.0, 3: 13.2, 4: 14.2, 5: 15.2, 6: 16.2, 7: 17.2, 8: 18.2, 9: 19.2},
        'kras_1': {0: 0.0, 1: 110.0, 2: 120.0, 3: 130.0, 4: 140.0, 5: 150.0, 6: 160.0, 7: 170.0, 8: 180.0, 9: 190.0},
        'kras_2': {0: 0.1, 1: 110.1, 2: 120.1, 3: 130.1, 4: 140.1, 5: 150.1, 6: 160.1, 7: 170.1, 8: 180.1, 9: 190.1},
        'kras_3': {0: 0.2, 1: 110.2, 2: 120.2, 3: 130.2, 4: 140.2, 5: 150.2, 6: 160.2, 7: 170.2, 8: 180.2, 9: 190.2}
        }
    
    expectedMatrixPDF = pd.DataFrame( expectedMatrixDict )
    
    # ################################################################################
    # def setUp(self):
    #     pass
    #
    # ################################################################################
    # def tearDown(self):
    #     pass

    ################################################################################
    def testloadSalmonReadsTable(self):
        self.logger.info( "BEGIN" )
        
        self.logger.info("expectedMatrixPDF:\n{}".format(self.expectedMatrixPDF))
        rawCountsSparkDF = CountMatrix(self.spark).loadSalmonReadsTable( self.fileList, self.sampleNames)
        retPDF = rawCountsSparkDF.toPandas()
        pd.testing.assert_frame_equal( self.expectedMatrixPDF,  retPDF)
        
        #retPDF.to_csv("./data/numReadsMatrix.tsv", index=False, sep="\t")
        
        self.logger.info( "END\n" )
        
    ################################################################################
    def testloadSalmonReadsTableWithRowId(self):
        self.logger.info( "BEGIN" )
        
        rawCountsSparkDF = CountMatrix(self.spark).loadSalmonReadsTableWithRowId( self.fileList, self.sampleNames)
        rawCountsSparkDF.show()
        pd.testing.assert_frame_equal( self.expectedMatrixPDF, rawCountsSparkDF.toPandas() )
        
        self.logger.info( "END\n" )
        

################################################################################
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()