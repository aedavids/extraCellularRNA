#!/usr/local/bin/python3
# encoding: utf-8
'''
bigDataDeseq.parseSalmonReads -- shortdesc

bigDataDeseq.parseSalmonReads is a description

It defines classes_and_methods

@author:     Andrew E. Davidson

@copyright:  2021 organization_name. All rights reserved.

@license:    license

@contact:    aedavids@ucsc.edu
@deffield    updated: Updated
'''


from   argparse import ArgumentParser
from   argparse import RawDescriptionHelpFormatter
# from bigDataDeseq.estimateScalingFactors import EstimateScalingFactors
from bigDataDeseq.countMatrix import CountMatrix

import pandas as pd
from   pyspark.sql import SparkSession


__all__ = []
__version__ = 0.1
__date__ = '2021-11-08'
__updated__ = '2021-11-08'
    
###############################################################################
class CommandLine( object ):
    '''
    Handle the command line, usage and help requests.
    '''

    def __init__( self, inOpts=None ):
        '''
        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''

        program_version = "v%s" % __version__
        program_build_date = str( __updated__ )
        program_version_message = '%%(prog)s %s (%s)' % ( program_version, program_build_date )
        program_shortdesc = "creates a single count matrix file from a set of salmon quant files and calculates \n"\
        + "the DESeq2 equivalent scaling factors"
        # program_shortdesc = __import__( '__main__' ).__doc__.split( "\n" )[1]
        # print("WTF: {}".format(__import__( '__main__' ).__doc__.split( "\n" )))
        # print("AEDWIP program_shortdesc: " + program_shortdesc + "XXXX")

        program_license = '''%s

      Created by Andrew E. Davidson on %s.
      Copyright 2020 organization_name. All rights reserved.

      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0

      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.

    USAGE
    ''' % ( program_shortdesc, str( __date__ ) )

        self.parser = ArgumentParser( description=program_license, formatter_class=RawDescriptionHelpFormatter )
        self.parser.add_argument( '-v', '--version', action='version', version=program_version_message )

        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

        # metavar
        # see https://stackoverflow.com/questions/26626799/pythons-argument-parser-printing-the-argument-name-in-upper-case
        # self.requiredArg.add_argument( '-m', '--mappingCSV', required=True, default=None, metavar="",
        #                                       action='store', help='mapping file in csv format with columns txId and geneId' )
        self.requiredArg.add_argument( '-q', '--quantFilesCSV', required=True, default=None, metavar="",
                                              action='store', help='file in csv format with two columns sampleName, and source. The source can be unix file path or url' )

        self.requiredArg.add_argument( '-o', '--outputDir', required=True, default=None, metavar="",
                                             action='store', help='path to unix directory or URL to GCP bucket, or AWS s3' )

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args( inOpts )

########################################################################
def main( inComandLineArgsList=None ):
    '''
    by default command line arguments will be parsed from sys.argv.
    use inComandLineArgsList to set arguments from unit tests or jupyter notebooks
    '''
    
    # GCP dataproc you have to set spark.driver.memory when you create the cluster
    # conf = SparkConf()
    # conf.set('spark.driver.memory', '16g')
     
    spark = SparkSession\
                .builder\
                .appName("countMatrixCLI")\
                .getOrCreate()
                #.config("spark.driver.memory", "15g")     .getOrCreate()   
        
    #
    # https://medium.com/@lubna_22592/building-production-pyspark-jobs-5480d03fd71e
    # initialize  logger for yarn cluster logs
    #
    log4jLogger = spark.sparkContext._jvm.org.apache.log4j
    logger = log4jLogger.LogManager.getLogger(__name__)
    logger.info("pyspark script logger initialized")   
    
    #
    # debug dump configuration
    #
    sparkConfig = spark.sparkContext.getConf().getAll()
    for item in sparkConfig: 
        logger.warn("sparkConfig: {}".format(item))
        print("sparkConfig: {}".format(item))
            
    #print("AEDWIP goes to stdout. if deployment mode = 'client' i.e. will stream to job detail page")
        
    if inComandLineArgsList is None:
        # parse arguments from sys.argv
        cli = CommandLine()
    else:
        cli = CommandLine( inComandLineArgsList )
            
    # python logger uses logger.warning()
    # spark uses log4j. logger.warning() generates reflection error
    logger.warn("arguments:\n {}".format(cli.args))
    
    # txId2GeneIdFile = cli.args.mappingCSV
    
    # get the sample names 
    quantFiles = cli.args.quantFilesCSV
    quantFilesPDF = pd.read_csv( quantFiles, index_col=False )

    logger.debug("col names:{}".format(quantFilesPDF.columns))
    logger.debug("quantFilesPDF:\n{}".format(quantFilesPDF))
    
    sampleNameList = list( quantFilesPDF.loc[:, "sampleName"] ) 
    # we can use the number of samples as sanity check
    logger.warn("number of samples to process:{}".format(len(sampleNameList)))
    
    # get file list
    fileList = list( quantFilesPDF.loc[:,"source"] )
    logger.debug("fileList:\n{}".format(fileList))

    outputDir = cli.args.outputDir
    
    # run
    logger.info("main() start execution")
    cm = CountMatrix(spark, log4jLogger)
    
    # countMatrixSparkDF = cm.loadSalmonReadsTableWithRowId(fileList, sampleNameList)
    countMatrixSparkDF = cm.loadSalmonReadsTable(fileList, sampleNameList)
    # AEDWIP TODO use CountMatrix, constructor args have changed 
    # AEDWIP change the name of this file it loads the counts and estimates
    # esf = EstimateScalingFactors( spark, fileList, sampleNameList, txId2GeneIdFile, log4jLogger )
    # retScalingFactorsDF, retCountDF = esf.run()
    
    # save
    if outputDir[-1] == "/" :
        outputDir = outputDir[:-1]
         
    # outputDir = outputDir + "/" + "preprocessData"
    # outFileESF = outputDir + "/" + "estimatedScalingFactors"
    outfileCount = outputDir + "/" + "countMatrixCLI.out/counts.tsv"   

    # spark uses log4j. logger.warning() generates reflection err0r
    # logger.warn("writing scalingFactors to {}".format(outFileESF ) )
    # retScalingFactorsDF.coalesce(1).write.csv( outFileESF, mode='overwrite',  header=True)       
    
    # spark uses log4j. logger.warning() generates reflection err0r    
    logger.warn("writing count matrix to {}".format(outfileCount ) ) 
    # do not coalesce. Spark can read/write parts in parallel
    # write a separate jobs to batch, combine, ... for Terra/DESeq    
    
    
    # coalesce(1) creates a singe part file
    countMatrixSparkDF.coalesce(1).write.csv( outfileCount, 
                                              mode='overwrite',
                                              sep="\t",
                                               header=True)

                 
########################################################################
if __name__ == "__main__":
    '''
    aedwip 0
    aedwip short discription
    '''
    
    # import sys;sys.argv = ['', 'Test.testName']
    main()

