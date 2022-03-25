'''
Created on Dec 22, 2021

@author: andrewdavidson
'''
# use findspark to inject pyspark into process
import findspark
import os
from absl.testing.absltest import _makedirs_exist_ok
# import bigDataDeseq

SPARK_HOME_ENV = os.getenv( "SPARK_HOME" )
if SPARK_HOME_ENV:
    # logger.info("using SPARK_HOME environment variable:{}".format(SPARK_HOME))
    findspark.init()
else:
    SPARK_HOME_HACK = "../../spark-3.1.2-bin-hadoop3.2"
    findspark.init( SPARK_HOME_HACK )

from   bigDataDeseq.parseSalmonReads import ParseSalmonReads
import logging
# import numpy as np
import pandas as pd
import pathlib as pl
from   pyspark.sql import SparkSession
from   bigDataDeseq.setupLogging import setupLogging
import shutil
import unittest


################################################################################
class TestParseSalmonReads( unittest.TestCase ):
    '''
    test cases are not independent. The will be run in alphabetic order
    '''

    # use class variable to make intermediate results avaliable to
    # down stream tests.
    globalDict = dict()

    spark = SparkSession\
                .builder\
                .appName( "TestEstimatedScalingFactors" )\
                .getOrCreate()
                    # .config("spark.driver.memory", "15g")     .getOrCreate()

    fileList = [
                "data/ctrl_1.quant.sf",
                "data/ctrl_2.quant.sf",
                "data/ctrl_3.quant.sf",
                "data/kras_1.quant.sf",
                "data/kras_2.quant.sf",
                "data/kras_3.quant.sf"
        ]

        # test we can work with sample names like GTEX-1117F-0426-SM-5EGHI
    sampleNames = [ "ctrl-1", "ctrl_2", "ctrl_3", "kras_1", "kras_2", "kras_3" ]

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

    tmp = pl.Path( "./tmp" )
    # ugly hack
    tmp.mkdir( exist_ok=True )
    logger.info( "pwd:{}".format( os.getcwd() ) )
    logger.info( "tmp:{}".format( tmp ) )
    shutil.rmtree( tmp )

    ################################################################################
    def setUp( self ):
        pass

    ################################################################################
    def tearDown( self ):
        pass

    ################################################################################
    def test_a_ParseSalmonReads( self ):
        self.logger.info( "BEGIN" )

        psr = ParseSalmonReads( self.spark, self.fileList, self.sampleNames )
        self.globalDict['psr'] = psr
        # resultsList list of dataframes. One for each quant file. contains 2 cols
        #  'name' and "numreads'
        resultsList = psr.selectCounts()
        self.globalDict['selectCounts'] = resultsList

        # tmp = pl.Path("./tmp/selectCounts")
        testTmp = self.tmp.joinpath( "selectCounts" )
        testTmp.mkdir( parents=True, exist_ok=True )

        for i in range( len( resultsList ) ):
            df = resultsList[i]
            sampleName = self.sampleNames[i]
            outfile = testTmp.joinpath( sampleName )
            self.logger.info( "outfile:{}".format( outfile.as_posix() ) )
            df.write.csv( outfile.as_posix(), mode='overwrite', header=True )

        # TODO add assert
        self.logger.info( "END" )

    ################################################################################
    def test_b_split( self ):
        self.logger.info( "BEGIN" )

        psr = self.globalDict['psr']
        selectCountsDFs = self.globalDict['selectCounts']
        # listOfPartDFLists = [None] * len(selectCountsDFs)
        # self.globalDict['listOfPartDFLists'] = listOfPartDFLists

        # we expect each df to be split with 3, 3, 3, and 1 data rows
        numParts = 3

        # tmp = pl.Path("./tmp/split")
        testTmp = self.tmp.joinpath( "split" )
        testTmp.mkdir( parents=True, exist_ok=True )

        # bash
        # for i in `seq 0 3`; do echo "\n $i"; cut -d , -f 1 tmp/split/${i}/ctrl_2/part*; done

        for i in range( len( selectCountsDFs ) ):
            df = selectCountsDFs[i]
            # we expect each df to be split with 3, 3, and 1 data rows
            resultsList = psr.split( df, numParts )
            # listOfPartDFLists[i] = resultsList
            sampleName = self.sampleNames[i]

            for j in range( len( resultsList ) ):
                partDF = resultsList[j]
                self.logger.info( "j:{} partDF".format( j ) )
                # partDF.show()
                outdir = testTmp.joinpath( str( j ) )
                outdir.mkdir( parents=True, exist_ok=True )
                outfile = outdir.joinpath( sampleName )
                self.logger.info( "outfile:{}".format( outfile.as_posix() ) )
                partDF.write.csv( outfile.as_posix(), mode='overwrite', header=True )

        self.logger.info( "END\n" )

    ################################################################################
    def test_c_joinParts( self ):
        self.logger.info( "BEGIN" )

        # find all the part files
        # test spark driver load from cloud storage
        splitDir = self.tmp.joinpath( "split" )
        partFileDict = dict()
        for partDir in splitDir.iterdir():
            if partDir.is_file():
                self.logger.info( "is file:{}".format( partDir ) )
            elif partDir.is_dir():
                self.logger.info( "is dir:{}".format( partDir ) )
                key = str( partDir ).split( "/" )[-1]
                self.logger.info( "key:{}".format( key ) )
                # partFileDict[key] = list()
                for sampleDir in partDir.iterdir():
                    if sampleDir.is_file():
                        self.logger.info( "sampleDir is file:{}".format( sampleDir ) )
                    elif sampleDir.is_dir():
                        sampleName = str( sampleDir ).split( "/" )[-1]
                        self.logger.info( "sampleName:{}".format( sampleName ) )
                        self.logger.info( "sampleDir is dir:{}".format( sampleDir ) )
                        for sampleFile in sampleDir.iterdir():
                            sampleFileStr = str( sampleFile )
                            # self.logger.info("aedwip sampleFileStr:{}".format(sampleFileStr))
                            if sampleFile.is_file() and sampleFileStr.endswith( ".csv" ):
                                self.logger.info( "sampleFile is file:{}".format( sampleFileStr ) )
                                if not key in partFileDict:
                                    partFileDict[key] = list()
                                t = ( sampleName, sampleFileStr )
                                partFileDict[key].append( t )

        self.logger.info( "partFileDict:\n{}".format( partFileDict ) )
        # {'0': [
        #    ('ctrl_2', 'tmp/split/0/ctrl_2/part-00000-59fdd373-dae0-486a-9c57-5bf992ca17df-c000.csv'),
        #    ('ctrl_3', 'tmp/split/0/ctrl_3/part-00000-432339d0-5754-40ef-b515-c366b9d7f397-c000.csv'),
        #    ('kras_1', 'tmp/split/0/kras_1/part-00000-000e9299-0f85-4102-b08a-b69a1dffa461-c000.csv'),
        #    ('kras_2', 'tmp/split/0/kras_2/part-00000-372b8310-e89f-4332-a946-a71c50d9782a-c000.csv'),
        #    ('kras_3', 'tmp/split/0/kras_3/part-00000-3dc6e996-8a7b-4019-878c-21a617e8806c-c000.csv'),
        #    ('ctrl-1', 'tmp/split/0/ctrl-1/part-00000-2c5fc27b-94d2-4a06-92ba-eb7a97fe3f77-c000.csv')],
        #
        # '1': [
        #   ('ctrl_2', 'tmp/split/1/ctrl_2/part-00000-c7c2f9f1-d0e3-491d-8655-d3af10682331-c000.csv'),
        #   ('ctrl_3', 'tmp/split/1/ctrl_3/part-00000-46cd4e42-d1a8-413a-8506-5f9cd4baf3d0-c000.csv'),
        #   ('kras_1', 'tmp/split/1/kras_1/part-00000-0146b231-1cb4-45c2-b945-3f3cb1034cf1-c000.csv'),
        #   ('kras_2', 'tmp/split/1/kras_2/part-00000-21b37bf6-627c-4792-8a7b-1df803b44683-c000.csv'),
        #   ('kras_3', 'tmp/split/1/kras_3/part-00000-673537f9-77ca-4aa8-9417-383a612ee7bf-c000.csv'),
        #   ('ctrl-1', 'tmp/split/1/ctrl-1/part-00000-f781c11f-c347-4491-a50e-dc16c2f5a237-c000.csv')],
        # ...

        expectedDict0 = {'Name': {0: 'txId_1', 1: 'txId_10', 2: 'txId_2'},
                        'ctrl-1': {0: 0.0, 1: 19.0, 2: 11.0},
                        'ctrl_2': {0: 0.1, 1: 19.1, 2: 11.1},
                        'ctrl_3': {0: 0.2, 1: 19.2, 2: 11.2},
                        'kras_1': {0: 0.0, 1: 190.0, 2: 110.0},
                        'kras_2': {0: 0.1, 1: 190.1, 2: 110.1},
                        'kras_3': {0: 0.2, 1: 190.2, 2: 110.2}}

        expectedDict1 = {'Name': {0: 'txId_3', 1: 'txId_4', 2: 'txId_5'},
                        'ctrl-1': {0: 12.0, 1: 13.0, 2: 14.0},
                        'ctrl_2': {0: 12.1, 1: 13.1, 2: 14.1},
                        'ctrl_3': {0: 0.0, 1: 13.2, 2: 14.2},
                        'kras_1': {0: 120.0, 1: 130.0, 2: 140.0},
                        'kras_2': {0: 120.1, 1: 130.1, 2: 140.1},
                        'kras_3': {0: 120.2, 1: 130.2, 2: 140.2}}

        expectedDict2 = {'Name': {0: 'txId_6', 1: 'txId_7', 2: 'txId_8'},
                        'ctrl-1': {0: 15.0, 1: 16.0, 2: 17.0},
                        'ctrl_2': {0: 15.1, 1: 16.1, 2: 17.1},
                        'ctrl_3': {0: 15.2, 1: 16.2, 2: 17.2},
                        'kras_1': {0: 150.0, 1: 160.0, 2: 170.0},
                        'kras_2': {0: 150.1, 1: 160.1, 2: 170.1},
                        'kras_3': {0: 150.2, 1: 160.2, 2: 170.2}}

        expectedDict3 = {'Name': {0: 'txId_9'},
                        'ctrl-1': {0: 18.0},
                        'ctrl_2': {0: 18.1},
                        'ctrl_3': {0: 18.2},
                        'kras_1': {0: 180.0},
                        'kras_2': {0: 180.1},
                        'kras_3': {0: 180.2}}

        expectedResults = { 0: pd.DataFrame(expectedDict0), 
                            1: pd.DataFrame(expectedDict1), 
                            2: pd.DataFrame(expectedDict2), 
                            3: pd.DataFrame(expectedDict3) }

        testTmp = self.tmp.joinpath( "joinedParts" )
        testTmp.mkdir( parents=True, exist_ok=True )
        
        psr = self.globalDict['psr']
        for partId, listOfTups in partFileDict.items():
            partId = int(partId)
            dfList = [None] * len( listOfTups )
            for i in range( len( listOfTups ) ):
                t = listOfTups[i]
                sampleName, filePath = t
                partSchema = "`Name` STRING, `{}` DOUBLE".format( sampleName )
                df = self.spark.read.load( filePath,
                           format="csv",
                           schema=partSchema,
                           header="true" )
                dfList[i] = df

            partCountDF = psr.joinParts( dfList )
            self.logger.info( "partId:{} ".format( partId ) )
            partCountDF.show()
            
            #
            # test
            #
            partCountPDF = partCountDF.toPandas()
            expectedPDF = expectedResults[ partId ]
            self.logger.info("partId:{} expectedPDF:\n{}".format(partId, expectedPDF))
            
            pd.testing.assert_frame_equal( partCountPDF, expectedPDF, check_dtype=False ) 
            
            #
            # save
            #      
            # outdir = testTmp.joinpath( str( partId ) )
            # outdir.mkdir( parents=True, exist_ok=True )
            outfile = testTmp.joinpath( str(partId) )
            self.logger.info( "outfile:{}".format( outfile.as_posix() ) )
            partCountDF.write.csv( outfile.as_posix(), mode='overwrite', header=True )

            print( " " )

        self.logger.info( "END\n" )
        
    ################################################################################
    def test_d_unionParts( self ):
        self.logger.info( "BEGIN" )
        
        # test spark driver load from cloud storage
        jointPartDFList = []
        joinedPartsDir = self.tmp.joinpath( "joinedParts" )
        
        # create schema
        sortedSampleNames = sorted(self.sampleNames )
        schema = "`Name` STRING"
        for sampleName in sortedSampleNames:
            schema = "{}, `{}` DOUBLE".format(schema, sampleName)
            
        self.logger.info("schema:{}".format(schema))
        
        expetectPDF = pd.DataFrame(
            {'Name': {0: 'txId_1', 1: 'txId_10', 2: 'txId_2', 3: 'txId_3', 4: 'txId_4', 
                      5: 'txId_5', 6: 'txId_9', 7: 'txId_6', 8: 'txId_7', 9: 'txId_8'}, 
             'ctrl-1': {0: 0.0, 1: 19.0, 2: 11.0, 3: 12.0, 4: 13.0, 5: 14.0, 6: 18.0, 
                        7: 15.0, 8: 16.0, 9: 17.0}, 
             'ctrl_2': {0: 0.1, 1: 19.1, 2: 11.1, 3: 12.1, 4: 13.1, 5: 14.1, 6: 18.1, 
                        7: 15.1, 8: 16.1, 9: 17.1}, 
             'ctrl_3': {0: 0.2, 1: 19.2, 2: 11.2, 3: 0.0, 4: 13.2, 5: 14.2, 6: 18.2, 
                        7: 15.2, 8: 16.2, 9: 17.2}, 
             'kras_1': {0: 0.0, 1: 190.0, 2: 110.0, 3: 120.0, 4: 130.0, 5: 140.0, 6: 180.0, 
                        7: 150.0, 8: 160.0, 9: 170.0}, 
             'kras_2': {0: 0.1, 1: 190.1, 2: 110.1, 3: 120.1, 4: 130.1, 5: 140.1, 6: 180.1, 
                        7: 150.1, 8: 160.1, 9: 170.1}, 
             'kras_3': {0: 0.2, 1: 190.2, 2: 110.2, 3: 120.2, 4: 130.2, 5: 140.2, 6: 180.2, 
                        7: 150.2, 8: 160.2, 9: 170.2} }
            )
        
        # load joined parts data frames
        for partDir in joinedPartsDir.iterdir():
            if partDir.is_dir():
                df = self.spark.read.load( partDir.as_posix( ),
                           format="csv",
                           schema=schema,
                           header="true" )
                jointPartDFList.append( df )
        
        # combine them into a single data frame
        psr = self.globalDict['psr']
        resultDF = psr.unionParts(jointPartDFList)
        self.logger.info("resultsDF")
        resultDF.show()
        resultPDF = resultDF.toPandas()
        # print("\ndict")
        # print(resultPDF.to_dict())
        
        #
        # test
        #
        pd.testing.assert_frame_equal( resultPDF, expetectPDF, check_dtype=False ) 

        self.logger.info( "END\n" )

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
