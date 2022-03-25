#!/usr/local/bin/python3
# encoding: utf-8

'''
Created on Dec 22, 2021

@author: andrewdavidson
'''

from   argparse import ArgumentParser
from   argparse import RawDescriptionHelpFormatter
from   bigDataDeseq.parseSalmonReads import ParseSalmonReads
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
        program_shortdesc = "Selects the 'name' and 'numreads' column from a salmon quant file and saves to output \n"
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
        self.requiredArg.add_argument( '-n', '--numParts', required=True, default=None, metavar="", type=int,
                                              action='store', help='number data partitions to be created' )

        self.requiredArg.add_argument( '-q', '--quantFilesCSV', required=True, default=None, metavar="",
                                              action='store', help='file in csv format with two columns sampleName, and source. The source can be unix file path or url. should be output from output by parseSalmonReadsSelectCountsCLI' )

        self.requiredArg.add_argument( '-o', '--outputDir', required=True, default=None, metavar="",
                                             action='store', help='path to unix directory or URL to GCP bucket, or AWS s3' )
        
        self.parser.add_argument( '-t', '--tag', default=" ", metavar="",
                                        action='store', help='add tag to log message' )

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args( inOpts )
            
################################################################################
def main(inComandLineArgsList=None ):
    '''
    by default command line arguments will be parsed from sys.argv.
    use inComandLineArgsList to set arguments from unit tests or jupyter notebooks
    '''
     
    spark = SparkSession\
                .builder\
                .appName("parseSalmonReads")\
                .getOrCreate()
                #.config("spark.driver.memory", "15g")     .getOrCreate()
                
    #
    # https://medium.com/@lubna_22592/building-production-pyspark-jobs-5480d03fd71e
    # initialize  logger for yarn cluster logs
    #
    log4jLogger = spark.sparkContext._jvm.org.apache.log4j
    logger = log4jLogger.LogManager.getLogger(__name__)
    logger.info("pyspark script logger initialized")  
    
    if inComandLineArgsList is None:
        # parse arguments from sys.argv
        cli = CommandLine()
    else:
        cli = CommandLine( inComandLineArgsList )
            
    # python logger uses logger.warning()
    # spark uses log4j. logger.warning() generates reflection error
    logger.warn("arguments:\n {}".format(cli.args))

    logTag = cli.args.tag

    # get the sample names 
    quantFiles = cli.args.quantFilesCSV
    quantFilesPDF = pd.read_csv( quantFiles, index_col=False )

    logger.warn("{} quantFilesPDF col names:{}".format(logTag, quantFilesPDF.columns))
    logger.debug("quantFilesPDF:\n{}".format(quantFilesPDF))
    
    sampleNameList = list( quantFilesPDF.loc[:, "sampleName"] ) 
    # we can use the number of samples as sanity check
    logger.warn("{} number of samples to process:{}".format(logTag, len(sampleNameList)))
    
    # get file list
    fileList = list( quantFilesPDF.loc[:,"source"] )    
    
    numParts = cli.args.numParts
    
    
    # get the output directory
    outputDir = cli.args.outputDir
    if outputDir[-1] == "/" :
        outputDir = outputDir[-1]
            
    # make it easy to find output
    outputDir = outputDir + "/parseSalmonReadsSplitCLI.out"

    # run
    psr = ParseSalmonReads( spark, None, None, log4jLogger ) # TODO fix ugly hack
    
    for i in range( len(fileList) ):
        filePath = fileList[i] 
        sampleName = sampleNameList[i]
        
        if filePath[-1] == "/" :
            t = filePath[:-1]
            #logger.warn("t:{}".format(t))
            filePath = t
            
        # expected file format
        # Name,ctrl-1
        # txId_1,0.0
        # txId_2,11.0
        # txId_3,12.0
        
        schema = "`Name` STRING, `{}` DOUBLE".format(sampleName)
        logger.warn("AEDWIP load filePath:{}".format(filePath))
        df = spark.read.load( filePath,
                           format="csv",
                           schema=schema,
                           header="true" )
        
        # split returns a list of dataframes
        # each data frame is of the form
        # Name,ctrl-1
        # txId_1,0.0
        # txId_10,19.0
        # txId_2,11.0
        # resultsList = psr.split(df, numberOfSplits=numParts)
        psr.split(df, numParts, sampleName, outputDir, logTag)

        
        # save
        # logger.warn("AEDWIP sampleName:{} len():{} ".format(sampleName, len( resultsList )))
        # for j in range( len( resultsList ) ):
        #     df = resultsList[j]
        #     outfile = outputDir + "/" + str(j) + "/" + sampleName 
        #     msg = "{} outfile:{}".format( logTag, outfile )
        #     logger.warn( msg )
        #     df.write.csv( outfile, mode='overwrite', header=True )

################################################################################
if __name__ == '__main__':
    main()