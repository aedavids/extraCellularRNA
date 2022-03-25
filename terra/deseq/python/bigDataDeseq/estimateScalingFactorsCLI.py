#!/usr/local/bin/python3
# encoding: utf-8
'''
bigDataDeseq.estimateScalingFactorsCLI -- Creates a Count Matrix and estimated scaling factors that can be used with DESeq2

bigDataDeseq.estimateScalingFactorsCLI is a command line wrapper around EstimatedScalingFactors

It defines classes_and_methods

@author:     Andrew E. Davidson

@copyright:  2022 organization_name. All rights reserved.

@license:    license

@contact:    aedavids@ucsc.edu
@deffield    updated: Updated
'''

import sys
import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from bigDataDeseq.estimateScalingFactors import EstimateScalingFactors

from   pyspark.sql import SparkSession


__all__ = []
__version__ = 0.1
__date__ = '2022-01-19'
__updated__ = '2022-01-19'

DEBUG = 0
TESTRUN = 0
PROFILE = 0
    
###############################################################################
class CommandLine( object ):
    '''
    Handle the command line, usage and help requests.
    '''

    def __init__( self, inOpts=None ):
        self.program_name = os.path.basename(sys.argv[0])
        program_version = "v%s" % __version__
        program_build_date = str(__updated__)
        program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
        program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
        program_license = '''%s
    
          Created by Andrew E. Davidson on %s.
          Copyright 2022 organization_name. All rights reserved.
        
          Licensed under the Apache License 2.0
          http://www.apache.org/licenses/LICENSE-2.0
        
          Distributed on an "AS IS" basis without warranties
          or conditions of any kind, either express or implied.
        
        USAGE
        ''' % (program_shortdesc, str(__date__))

        # Setup argument parser
        self.parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )
        self.requiredArg.add_argument( '-m', '--mappingCSV', required=True, default=None, metavar="",
                                              action='store', help=("a unix file path or URL to GCP bucket, or AWS s3"
                                                                    " to a mapping file in csv format with columns txId and geneId"
                                                                    ) )
        
        self.requiredArg.add_argument( '-r', '--readsTSV', required=True, default=None, metavar="",
                                              action='store', help=("a unix file path or URL to GCP bucket, or AWS s3"
                                                                    " to a tab separated file. The first column is 'Name' this is the"
                                                                    " 'Name' values in your salmon quant.sf file" 
                                                                    " each of the remaining columns is the numReads column"
                                                                    " from the quant.sf files. The column name should be"
                                                                    " the sample name"
                                                                    ))
                                                   
        
        self.requiredArg.add_argument( '-o', '--outputDir', required=True, default=None, metavar="",
                                             action='store', help='path to unix directory or URL to GCP bucket, or AWS s3' )
        
        # optional arguments
        self.parser.add_argument('-v', '--version', action='version', version=program_version_message)        

        #
        # Process arguments
        #
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args( inOpts )
        
########################################################################
def main(inComandLineArgsList=None ): 
    '''
    by default command line arguments will be parsed from sys.argv.
    use inComandLineArgsList to set arguments from unit tests or jupyter notebooks
    '''
    try:
        # init spark at very beginning so we can get access to log4j
        spark = SparkSession\
                    .builder\
                    .appName("estimatedScalingFactors")\
                    .getOrCreate()    
                    
        #
        # https://medium.com/@lubna_22592/building-production-pyspark-jobs-5480d03fd71e
        # initialize  logger for yarn cluster logs
        #
        log4jLogger = spark.sparkContext._jvm.org.apache.log4j
        logger = log4jLogger.LogManager.getLogger(__name__)
        logger.warn("pyspark script logger initialized")   
                     
        # parse cli args
        if inComandLineArgsList is None:
            # parse arguments from sys.argv
            cli = CommandLine()
        else:
            cli = CommandLine( inComandLineArgsList )    
     
        # python logger uses logger.warning()
        # spark uses log4j. logger.warning() generates reflection error
        logger.warn("arguments:\n {}".format(cli.args)) 
        
        #
        # debug dump configuration only after parsing cli. 
        # we do not want to dump the spark config unless all required aruments 
        # are present
        #
        sparkConfig = spark.sparkContext.getConf().getAll()
        for item in sparkConfig: 
            logger.warn("sparkConfig: {}".format(item))
                  
        outputDirArg = cli.args.outputDir
        if outputDirArg[-1] == "/" :
            outputDirArg = outputDirArg[:-1]
        outputDir  = outputDirArg + "/" + cli.program_name     
               
        readsTSVFile = cli.args.readsTSV
        txId2GeneIdFile = cli.args.mappingCSV
        
        
        #
        # enable spark checkpoint 
        #
        sparkContext = spark.sparkContext
        #path must be accessible from the executor, check pointing is done by executor
        # hdfs://sparkCheckpoints" did not work
        #""
        #sparkContext.setCheckpointDir("hdfs://sparkCheckPointDir")                    
        #sparkContext.setCheckpointDir("gs://anvil-gtex-v8-hg38-edu-ucsc-kim-lab-spark/quant/sparkGTEx-train.out/sparkCheckPoints")
        checkpointDir = outputDir + "/checkpoinDir/"
        sparkContext.setCheckpointDir(checkpointDir)
        logger.warn("spark checkpoint dir:{}".format(checkpointDir))
            
        #
        # run
        #
        
        # it would be faster if we could provide the schema, how ever we do not know the col names
        readsSparkDF = spark.read.load( readsTSVFile, format="csv", sep="\t", header=True, inferSchema=True )
        
        # 1698 partitions created by default on cluster with 96 core and 2.8 Tb of memory 
        logger.warn("readsSparkDF.rdd.getNumPartitions():{}".format(readsSparkDF.rdd.getNumPartitions()))
        esf = EstimateScalingFactors( spark, txId2GeneIdFile )
        outFileGroupedByGene = outputDir + "/" + "countsGroupedByGene"         
        retScalingFactorsDF, retCountDF = esf.run(readsSparkDF, columnBatchSize=500, outFileGroupedByGene=outFileGroupedByGene)
        
        # save
        outFileESF = outputDir + "/" + "estimatedScalingFactors"
        outfileCount = outputDir + "/" + "counts"   
        
        # spark uses log4j. logger.warning() generates reflection err0r
        logger.warn("writing scalingFactors to {}".format(outFileESF ) )
        # use coalesce to create a single part file
        retScalingFactorsDF.coalesce(1).write.csv( outFileESF, mode='overwrite',  header=True)       
        
        # spark uses log4j. logger.warning() generates reflection err0r    
        logger.warn("writing count matrix to {}".format(outfileCount ) ) 
        # use coalesce to create a single part file
        retCountDF.coalesce(1).write.csv( outfileCount, mode='overwrite', header=True)
            
        return 0
    
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    
    except Exception as e:
        logger.error("ERROR main() {}".format(e))
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(cli.program_name) * " "
        sys.stderr.write(cli.program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help\n")
        return 2

########################################################################
if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'bigDataDeseq.estimateScalingFactorsCLI_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
        
    sys.exit(main())