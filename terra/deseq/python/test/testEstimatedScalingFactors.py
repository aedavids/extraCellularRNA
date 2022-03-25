'''
Created on Nov 9, 2021

@author: andrewdavidson
'''
# use findspark to inject pyspark into process
import findspark
import os
from bigDataDeseq.countMatrix import CountMatrix

SPARK_HOME_ENV = os.getenv( "SPARK_HOME" )
if SPARK_HOME_ENV:
    # logger.info("using SPARK_HOME environment variable:{}".format(SPARK_HOME))
    findspark.init()
else:
    SPARK_HOME_HACK = "../../spark-3.1.2-bin-hadoop3.2"
    findspark.init( SPARK_HOME_HACK )

import logging
import numpy as np
import pandas as pd
from   pyspark.sql import SparkSession
from   pyspark import SparkContext
import shutil

from   bigDataDeseq.setupLogging import setupLogging
import unittest

from bigDataDeseq.estimateScalingFactors import EstimateScalingFactors


################################################################################
class TestEstimatedScalingFactors( unittest.TestCase ):
    '''
    test cases are not independent. The will be run in alphabetic order
    '''
    spark = SparkSession\
                .builder\
                .appName( "TestEstimatedScalingFactors" )\
                .getOrCreate()
                    # .config("spark.driver.memory", "15g")     .getOrCreate()
                    
    sparkContext = spark.sparkContext
    #path must be accessible from the executor, check pointing is done by executor
    checkPointDir = "./sparkCheckPoints"
    sparkContext.setCheckpointDir(checkPointDir);                    

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

    # use class variable to make intermediate results available to
    # down stream tests.
    globalDict = dict()

    pd.set_option( 'display.max_rows', None )
    pd.set_option( 'display.max_columns', None )
    pd.set_option( 'display.width', None )
    pd.set_option( 'display.max_colwidth', -1 )

    # ################################################################################
    # def setUp(self):
    #     self.logger.info("BEGIN set up before every test")
    #
    #
    #     self.logger.info("END set up before every test")
    #     pass
    #
    # ################################################################################
    # def tearDown(self):
    #     self.logger.info("BEGIN clean up after every test")
    #
    #     self.logger.info("END   clean up after every test")
    #     pass
    
    @classmethod
    def tearDownClass(cls):
        cls.logger.warn("remove spark")
        shutil.rmtree( cls.checkPointDir )

    ################################################################################
    def test_A_groupByGeneAndSum( self ):
        self.logger.info( "BEGIN" )

        fileList = [
                    "data/ctrl_1.quant.sf",
                    "data/ctrl_2.quant.sf",
                    "data/ctrl_3.quant.sf",
                    "data/kras_1.quant.sf",
                    "data/kras_2.quant.sf",
                    "data/kras_3.quant.sf"
            ]

        # sampleNames = [ "ctrl_1", "ctrl_2", "ctrl_3", "kras_1", "kras_2", "kras_3" ]

        # test we can work with sample names like GTEX-1117F-0426-SM-5EGHI
        sampleNames = [ "ctrl-1", "ctrl_2", "ctrl_3", "kras_1", "kras_2", "kras_3" ]
        # txId2GeneIdFile = "data/txId2GeneId.csv"
        txId2GeneIdFile = "data/mockTxId2GeneId.csv"

        # load the numReads from the salmon quant.sf files
        esf = EstimateScalingFactors( self.spark, txId2GeneIdFile )
        self.globalDict['esf'] = esf

        rawCountsSparkDF = CountMatrix( self.spark ).loadSalmonReadsTable( fileList, sampleNames )
        rawCountsSparkDF.createOrReplaceTempView( "rawCounts" )
        self.logger.info( "\nrawCountsSparkDF.show()" )
        rawCountsSparkDF.show( truncate=False )

        expectedRawCountsDict = {
            'Name': {0: 'txId_1', 1: 'txId_2', 2: 'txId_3', 3: 'txId_4', 4: 'txId_5', 5: 'txId_6', 6: 'txId_7', 7: 'txId_8', 8: 'txId_9', 9: 'txId_10'},
            'ctrl-1': {0: 0.0, 1: 11.0, 2: 12.0, 3: 13.0, 4: 14.0, 5: 15.0, 6: 16.0, 7: 17.0, 8: 18.0, 9: 19.0},
            'ctrl_2': {0: 0.1, 1: 11.1, 2: 12.1, 3: 13.1, 4: 14.1, 5: 15.1, 6: 16.1, 7: 17.1, 8: 18.1, 9: 19.1},
            'ctrl_3': {0: 0.2, 1: 11.2, 2: 0.0, 3: 13.2, 4: 14.2, 5: 15.2, 6: 16.2, 7: 17.2, 8: 18.2, 9: 19.2},
            'kras_1': {0: 0.0, 1: 110.0, 2: 120.0, 3: 130.0, 4: 140.0, 5: 150.0, 6: 160.0, 7: 170.0, 8: 180.0, 9: 190.0},
            'kras_2': {0: 0.1, 1: 110.1, 2: 120.1, 3: 130.1, 4: 140.1, 5: 150.1, 6: 160.1, 7: 170.1, 8: 180.1, 9: 190.1},
            'kras_3': {0: 0.2, 1: 110.2, 2: 120.2, 3: 130.2, 4: 140.2, 5: 150.2, 6: 160.2, 7: 170.2, 8: 180.2, 9: 190.2}
            }
        expectedRawCountsPDF = pd.DataFrame( expectedRawCountsDict )
        pd.testing.assert_frame_equal( expectedRawCountsPDF, rawCountsSparkDF.toPandas() )

        # rawCountPDF = rawCountsSparkDF.toPandas()
        # self.logger.info("rawCountPDF dict:\n{}".format(rawCountPDF.to_dict()))

        # save countSparkDF, it is needed for test_B_
        # groupByGeneAndSum changes the row order. this does not matter
        countsSparkDF = esf._groupByGeneAndSum( rawCountsSparkDF )
        self.logger.info( "countsSparkDF.toPandas \n{}".format( countsSparkDF.toPandas() ) )
        # self.logger.info("countsSparkDF.toPandas.to_dict() \n{}".format(countsSparkDF.toPandas().to_dict()))

        self.globalDict["countsSparkDF"] = countsSparkDF

        #
        # test the grouped counts are correct
        #
        expectedCountDict = {
            'geneId': {0: 'gene_2', 1: 'gene_7', 2: 'gene_6', 3: 'gene_4', 4: 'gene_1', 5: 'gene_5', 6: 'gene_3', 7: 'gene_8'},
            'sum(ctrl-1)': {0: 11.0, 1: 16.0, 2: 15.0, 3: 13.0, 4: 0.0, 5: 14.0, 6: 12.0, 7: 54.0},
            'sum(ctrl_2)': {0: 11.1, 1: 16.1, 2: 15.1, 3: 13.1, 4: 0.1, 5: 14.1, 6: 12.1, 7: 54.300000000000004},
            'sum(ctrl_3)': {0: 11.2, 1: 16.2, 2: 15.2, 3: 13.2, 4: 0.2, 5: 14.2, 6: 0.0, 7: 54.599999999999994},
            'sum(kras_1)': {0: 110.0, 1: 160.0, 2: 150.0, 3: 130.0, 4: 0.0, 5: 140.0, 6: 120.0, 7: 540.0},
            'sum(kras_2)': {0: 110.1, 1: 160.1, 2: 150.1, 3: 130.1, 4: 0.1, 5: 140.1, 6: 120.1, 7: 540.3},
            'sum(kras_3)': {0: 110.2, 1: 160.2, 2: 150.2, 3: 130.2, 4: 0.2, 5: 140.2, 6: 120.2, 7: 540.5999999999999}}

        expectedCountPDF = pd.DataFrame( expectedCountDict )
        pd.testing.assert_frame_equal( expectedCountPDF, countsSparkDF.toPandas() )

        #
        # test counts are correct
        #
        expectedRowStr = "Row(geneId='gene_8', sum(ctrl-1)=54.0, sum(ctrl_2)=54.300000000000004, sum(ctrl_3)=54.599999999999994, sum(kras_1)=540.0, sum(kras_2)=540.3, sum(kras_3)=540.5999999999999)"
        retRow = countsSparkDF.filter( countsSparkDF.geneId == "gene_8" ).take( 1 )[0]
        retRowStr = str( retRow )
        self.logger.info( "expectedRowStr:{}".format( expectedRowStr ) )
        self.logger.info( "     retRowStr:{}".format( retRowStr ) )
        self.assertEqual( expectedRowStr, retRowStr, "_groupByGeneAndSum() failed" )

        self.logger.info( "END\n" )

    ################################################################################
    def test_B_convertCountsToInts( self ):
        self.logger.info( "BEGIN" )

        esf = self.globalDict["esf"]
        countsSparkDF = self.globalDict["countsSparkDF"]

        intSparkDF = esf._convertToLong( countsSparkDF )
        self.globalDict["intSparkDF"] = intSparkDF

        expectedPDFDict = {'geneId': {0: 'gene_2', 1: 'gene_7', 2: 'gene_6', 3: 'gene_4', 4: 'gene_1', 5: 'gene_5', 6: 'gene_3', 7: 'gene_8'},
                       'sum(ctrl-1)': {0: 11, 1: 16, 2: 15, 3: 13, 4: 0, 5: 14, 6: 12, 7: 54},
                       'sum(ctrl_2)': {0: 11, 1: 16, 2: 15, 3: 13, 4: 0, 5: 14, 6: 12, 7: 54},
                       'sum(ctrl_3)': {0: 11, 1: 16, 2: 15, 3: 13, 4: 0, 5: 14, 6: 0, 7: 55},
                       'sum(kras_1)': {0: 110, 1: 160, 2: 150, 3: 130, 4: 0, 5: 140, 6: 120, 7: 540},
                       'sum(kras_2)': {0: 110, 1: 160, 2: 150, 3: 130, 4: 0, 5: 140, 6: 120, 7: 540},
                       'sum(kras_3)': {0: 110, 1: 160, 2: 150, 3: 130, 4: 0, 5: 140, 6: 120, 7: 541}}

        expectedPDF = pd.DataFrame( expectedPDFDict )
        retPDF = intSparkDF.toPandas()

        # spark df is int32, pandas from dict is int64
        self.logger.info( "expectedPDF:\n{}\n".format( expectedPDF ) )
        self.logger.info( "retPDF:\n{}".format( retPDF ) )
        pd.testing.assert_frame_equal( retPDF, expectedPDF, check_dtype=False )

        self.logger.info( "END\n" )

    ################################################################################
    def test_C_convertCountsToLog( self ):
        self.logger.info( "BEGIN" )

        esf = self.globalDict["esf"]
        intSparkDF = self.globalDict["intSparkDF"]
        # intSparkDF.show()

        columnNames = intSparkDF.columns[1:]
        logCountsSparkDF = esf._calculateLogs( intSparkDF, columnNames )
        # logCountsSparkDF.show()
        self.globalDict["logCountsSparkDF"] = logCountsSparkDF

        # test
        testPDF = logCountsSparkDF.select( ['geneId', 'log(ctrl-1)', 'log(kras_3)'] ).toPandas()
        # print(testPDF)
        # print(testPDF.to_dict())

        expectedPDFDict = {
         'geneId': {0: 'gene_2', 1: 'gene_7', 2: 'gene_6', 3: 'gene_4', 4: 'gene_1', 5: 'gene_5', 6: 'gene_3', 7: 'gene_8'},
         'log(ctrl-1)': {0: 2.3978952727983707, 1: 2.772588722239781, 2: 2.70805020110221, 3: 2.5649493574615367, 4: 'nan', 5: 2.6390573296152584, 6: 2.4849066497880004, 7: 3.9889840465642745},
         'log(kras_3)': {0: 4.700480365792417, 1: 5.075173815233827, 2: 5.0106352940962555, 3: 4.867534450455582, 4: 'nan', 5: 4.941642422609304, 6: 4.787491742782046, 7: 6.293419278846481}
        }

        # print(expectedPDF)
        expectedPDF = pd.DataFrame( expectedPDFDict )
        selectRows = expectedPDF.loc[:, "geneId"] == "gene_1"
        expectedPDF.loc[selectRows, ["log(ctrl-1)", "log(kras_3)"] ] = np.NaN  # pd.NaN
        # print(expectedPDF)

        # expectedPDF

        pd.testing.assert_frame_equal( testPDF, expectedPDF, check_dtype=False )

        self.logger.info( "END\n" )

    ################################################################################
    def test_D_filter( self ):
        self.logger.info( "BEGIN" )

        logCountsSparkDF = self.globalDict["logCountsSparkDF"]
        filteredDF = logCountsSparkDF.na.drop()
        self.globalDict["filteredDF"] = filteredDF

        # test
        expectedDict = {'geneId': {0: 'gene_2', 1: 'gene_7', 2: 'gene_6', 3: 'gene_4', 4: 'gene_5', 5: 'gene_8'}}

        expectedPDF = pd.DataFrame( expectedDict )
        filteredPDF = filteredDF.select( "geneId" ).toPandas()
        pd.testing.assert_frame_equal( filteredPDF , expectedPDF )

        self.logger.info( "END\n" )

    ################################################################################
    def test_E_rowSums( self ):
        self.logger.info( "BEGIN" )

        esf = self.globalDict["esf"]
        filteredDF = self.globalDict["filteredDF"]

        columns = filteredDF.columns[1:]
        rowSumsDF = esf.rowSums( filteredDF, columns, columnBatchSize=4 )
        self.globalDict["rowSumsDF"] = rowSumsDF

        # test
        pdf = rowSumsDF.select( ["geneId", "rowSum"] ).toPandas()
        # print(pdf.to_dict())

        expectedDict = {
            'geneId': {0: 'gene_2', 1: 'gene_7', 2: 'gene_6', 3: 'gene_4', 4: 'gene_5', 5: 'gene_8'},
            'rowSum': {0: 21.29512691577236, 1: 23.54328761242082, 2: 23.156056485595393, 3: 22.297451423751358, 4: 22.74209925667369, 5: 30.861858836324142}
        }

        expectedPDF = pd.DataFrame( expectedDict )

        pd.testing.assert_frame_equal( pdf , expectedPDF )

        self.logger.info( "END\n" )

    ################################################################################
    def test_F_rowMeans( self ):
        self.logger.info( "BEGIN" )

        rowSumsDF = self.globalDict["rowSumsDF"]

        n = len( rowSumsDF.columns ) - 2  # do not count geneId or rowSum columns
        rowMeansDF = rowSumsDF.withColumn( "rowMean", ( rowSumsDF.rowSum / n ) )
        self.globalDict["rowMeansDF"] = rowMeansDF

        # test
        pdf = rowSumsDF.select( ["geneId", "rowSum"] ).toPandas()
        # print(pdf.to_dict())

        expectedDict = {
            'geneId': {0: 'gene_2', 1: 'gene_7', 2: 'gene_6', 3: 'gene_4', 4: 'gene_5', 5: 'gene_8'},
            'rowSum': {0: 21.29512691577236, 1: 23.54328761242082, 2: 23.156056485595393, 3: 22.297451423751358, 4: 22.74209925667369, 5: 30.861858836324142}
        }

        expectedPDF = pd.DataFrame( expectedDict )

        pd.testing.assert_frame_equal( pdf , expectedPDF )

        self.logger.info( "END\n" )

    ################################################################################
    def test_G_subtractRowMeanFromLogCounts( self ):
        self.logger.info( "BEGIN" )

        esf = self.globalDict["esf"]
        rowMeansDF = self.globalDict["rowMeansDF"]

        columnNames = rowMeansDF.columns[1:-2]
        ratioDF = esf._subtractRowMeanFromLogCounts( rowMeansDF, columnNames )
        self.globalDict["ratioDF"] = ratioDF

        # test
        ratioPDF = ratioDF.toPandas()

        expectedDict = {
            'geneId': {0: 'gene_2', 1: 'gene_7', 2: 'gene_6', 3: 'gene_4', 4: 'gene_5', 5: 'gene_8'},
            'ctrl-1': {0:-1.151292546497023, 1:-1.151292546497022, 2:-1.151292546497022, 3:-1.151292546497023, 4:-1.1512925464970234, 5:-1.1546590928230822},
            'ctrl_2': {0:-1.151292546497023, 1:-1.151292546497022, 2:-1.151292546497022, 3:-1.151292546497023, 4:-1.1512925464970234, 5:-1.1546590928230822},
            'ctrl_3': {0:-1.151292546497023, 1:-1.151292546497022, 2:-1.151292546497022, 3:-1.151292546497023, 4:-1.1512925464970234, 5:-1.1363099541548856},
            'kras_1': {0: 1.151292546497023, 1: 1.1512925464970234, 2: 1.1512925464970234, 3: 1.1512925464970225, 4: 1.151292546497022, 5: 1.1479260001709637},
            'kras_2': {0: 1.151292546497023, 1: 1.1512925464970234, 2: 1.1512925464970234, 3: 1.1512925464970225, 4: 1.151292546497022, 5: 1.1479260001709637},
            'kras_3': {0: 1.151292546497023, 1: 1.1512925464970234, 2: 1.1512925464970234, 3: 1.1512925464970225, 4: 1.151292546497022, 5: 1.1497761394591244}}

        expectedPDF = pd.DataFrame( expectedDict )

        pd.testing.assert_frame_equal( ratioPDF , expectedPDF )

        self.logger.info( "END]\n" )

    ################################################################################
    def singleRowDFToNumpyArray( self, sparkDF ):
        pdf = sparkDF.toPandas()
        np = pdf.values[0]
        return np

    ################################################################################
    def test_H_median( self ):
        self.logger.info( "BEGIN" )

        # test even number of rows
        evenPDF = pd.DataFrame( {
            "a": [i * 1.0 for i in range( 6 )],
            "b": [i * 2.0 for i in range( 6 )],
        } )
        evenSparkDF = self.spark.createDataFrame( evenPDF )
        esf = self.globalDict["esf"]
        resultDF = esf.median( evenSparkDF, evenSparkDF.columns )
        resultMediansNP = self.singleRowDFToNumpyArray( resultDF )

        # test odd number of rows
        expectMedianNP = np.array( [2., 4.] )
        np.testing.assert_array_equal( resultMediansNP, expectMedianNP )

        oddPDF = pd.DataFrame( {
            "c": [i * 1.0 for i in range( 7 )],
            "d": [i * 2.0 for i in range( 7 )],
        } )
        oddSparkDF = self.spark.createDataFrame( oddPDF )
        resultDF = esf.median( oddSparkDF, oddSparkDF.columns )
        resultMediansNP = self.singleRowDFToNumpyArray( resultDF )
        expectMedianNP = np.array( [3., 6.] )
        np.testing.assert_array_equal( resultMediansNP, expectMedianNP )

        self.logger.info( "END\n" )

    ################################################################################
    def test_I_logScalingFactors( self ):
        self.logger.info( "BEGIN" )

        esf = self.globalDict["esf"]
        ratioDF = self.globalDict["ratioDF"]

        columnNames = ratioDF.columns[1:]
        logMedianDF = esf.median( ratioDF, columnNames )

        newColNames = [esf.getSampleNames( c ) for c in logMedianDF.columns]
        logScalingFactorsDF = logMedianDF.toDF( *newColNames )
        self.globalDict["logScalingFactorsDF"] = logScalingFactorsDF

        # test
        logScalingFactorsPDF = logScalingFactorsDF.toPandas()
        #     print(logScalingFactorsPDF.to_dict())

        expectedDict = {'ctrl-1': {0:-1.151292546497023}, 'ctrl_2': {0:-1.151292546497023},
                         'ctrl_3': {0:-1.151292546497023}, 'kras_1': {0: 1.1512925464970225},
                         'kras_2': {0: 1.1512925464970225}, 'kras_3': {0: 1.1512925464970225}}

        expectedPDF = pd.DataFrame( expectedDict )
        pd.testing.assert_frame_equal( logScalingFactorsPDF , expectedPDF )

        self.logger.info( "END\n" )

    ################################################################################
    def test_J_run( self ):
        self.logger.info( "BEGIN" )

        fileList = [
                    "data/ctrl_1.quant.sf",
                    "data/ctrl_2.quant.sf",
                    "data/ctrl_3.quant.sf",
                    "data/kras_1.quant.sf",
                    "data/kras_2.quant.sf",
                    "data/kras_3.quant.sf"
            ]

        sampleNames = [ "ctrl-1", "ctrl_2", "ctrl_3", "kras_1", "kras_2", "kras_3" ]
        txId2GeneIdFile = "data/mockTxId2GeneId.csv"

        # load the numReads from the salmon quant.sf files
        rawCountsSparkDF = CountMatrix( self.spark ).loadSalmonReadsTable( fileList, sampleNames )
        self.logger.info( "rawCountsSparkDF.show()" )
        rawCountsSparkDF.show()

        esf = EstimateScalingFactors( self.spark, txId2GeneIdFile )
        retScalingFactorsDF, retCountDF = esf.run( rawCountsSparkDF, columnBatchSize=4 )
        self.logger.info( "retScalingFactorsDF.show()" )
        retScalingFactorsDF.show( truncate=False )
        self.logger.info( "retCountDF" )
        retCountDF.show( truncate=False )

        # test scaling factors
        expectedScalingFactors = {
               'sampleName': {0: 'ctrl-1', 1: 'ctrl_2', 2: 'ctrl_3', 3: 'kras_1',
                              4: 'kras_2', 5: 'kras_3'},
            'scalingFactor': {0: 0.3162277660168379, 1: 0.3162277660168379, 2: 0.3162277660168379, 3: 3.162277660168378,
                              4: 3.162277660168378, 5: 3.162277660168378}}
        expectedScalingFactorPDF = pd.DataFrame( expectedScalingFactors )

        retScalingFactorsPDF = retScalingFactorsDF.toPandas()
        pd.testing.assert_frame_equal( retScalingFactorsPDF , expectedScalingFactorPDF )

        # test counts
        # version before we add row_number. difference is ordering of columns
        expectedCountsDict = {
            'geneId': {0: 'gene_2', 1: 'gene_7', 2: 'gene_6', 3: 'gene_4', 4: 'gene_1', 5: 'gene_5', 6: 'gene_3', 7: 'gene_8'},
            'ctrl-1': {0: 11, 1: 16, 2: 15, 3: 13, 4: 0, 5: 14, 6: 12, 7: 54},
            'ctrl_2': {0: 11, 1: 16, 2: 15, 3: 13, 4: 0, 5: 14, 6: 12, 7: 54},
            'ctrl_3': {0: 11, 1: 16, 2: 15, 3: 13, 4: 0, 5: 14, 6: 0, 7: 55},
            'kras_1': {0: 110, 1: 160, 2: 150, 3: 130, 4: 0, 5: 140, 6: 120, 7: 540},
            'kras_2': {0: 110, 1: 160, 2: 150, 3: 130, 4: 0, 5: 140, 6: 120, 7: 540},
            'kras_3': {0: 110, 1: 160, 2: 150, 3: 130, 4: 0, 5: 140, 6: 120, 7: 541}}

        # expectedCountsDict = {'geneId': {0: 'gene_1', 1: 'gene_2', 2: 'gene_3', 3: 'gene_4', 4: 'gene_5', 5: 'gene_6', 6: 'gene_7', 7: 'gene_8'},
        #                        'ctrl-1': {0: 0, 1: 11, 2: 12, 3: 13, 4: 14, 5: 15, 6: 16, 7: 54},
        #                        'ctrl_2': {0: 0, 1: 11, 2: 12, 3: 13, 4: 14, 5: 15, 6: 16, 7: 54},
        #                        'ctrl_3': {0: 0, 1: 11, 2: 0, 3: 13, 4: 14, 5: 15, 6: 16, 7: 55},
        #                        'kras_1': {0: 0, 1: 110, 2: 120, 3: 130, 4: 140, 5: 150, 6: 160, 7: 540},
        #                        'kras_2': {0: 0, 1: 110, 2: 120, 3: 130, 4: 140, 5: 150, 6: 160, 7: 540},
        #                        'kras_3': {0: 0, 1: 110, 2: 120, 3: 130, 4: 140, 5: 150, 6: 160, 7: 541}}

        expectedCountsPDF = pd.DataFrame( expectedCountsDict )
        self.logger.info( "expectedCountsPDF:\n{}".format( expectedCountsPDF ) )

        retCountPDF = retCountDF.toPandas()
        print( "AEDWIP !!!!!!!!!!!!!!!!!!" )
        print( retCountPDF.to_dict() )

        self.logger.info( "retCountPDF:\n{}".format( retCountPDF ) )
        pd.testing.assert_frame_equal( retCountPDF , expectedCountsPDF, check_dtype=False )

        self.logger.info( "END\n" )


################################################################################
if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
