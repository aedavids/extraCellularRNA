#
# findBiomarkers.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref:

# findBiomarkers display the doc string 
'''
findBiomarkers.py TODO doc string for cli
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from kimLabUtils.baseCommandLine import BBaseCommandLine
import logging
logger = logging.getLogger(__file__)
import os
import pandas as pd
import sys

__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2024-12-24'
__updated__ = '2024-12-24'

###############################################################################
class findBiomarkersCommandLine( BBaseCommandLine ):
    '''
    Handle the command line, usage and help requests.
    '''
    ###############################################################################
    def __init__( self, version, author, date, update ):
        '''
        Implement a parser to interpret the command line argv string using argparse.

        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''
        super().__init__(version, author, date, update )
        # self.author = author
        # self.program_version = version
        # self.program_build_date = str( update )
        # self.date = date
        #
        # self.program_version_message = '%%(prog)s %s (%s)' % ( self.program_version, self.program_build_date )
        # self.program_shortdesc = __import__( '__main__' ).__doc__.split( "\n" )[1]

    ###############################################################################
    def _build( self ):
        self.parser = ArgumentParser( description=self._getLicence(), formatter_class=RawDescriptionHelpFormatter )
    
        #
        # optional arguments
        #

        self.parser.add_argument( '-b', '--bioType', default=None, metavar="",
                                                      action='store', 
                                                      required=True,
                                                      help="select biomarkers from biotype. Default is to select from, all genes"
        )        
       
        self.parser.add_argument( '-c', '--countFilePath', default=None, 
                                                        metavar="",
                                                      action='store', 
                                                      required=False, # only required if biotype is specified
                                                      help="path to a csv file with gene counts and bio type. "
                                                        + "ex. '/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv'"
                                                        + "1st 3 column are gene_id, gene_name, gene_biotype there is a column for each sample"
        )
        

        # self.parser.add_argument( '-o', '--outDir', default=".", metavar="",
        #                                               action='store', 
        #                                               help="locaiton to write output files"
        # )

    
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

        # self.requiredArg.add_argument( '-o', '--outDir', default=".", metavar="",
        #                                               action='store', 
        #                                               help="locaiton to write output files mixture.tsv and signatureGenes.tsv"
        # )        

      
        # self.requiredArg.add_argument( '-m', '--metaDataFilePath', required=True, default=None, metavar="",
        #                                               action='store', 
        #                                               help="path to a csv file containing sample meta data in DESeq format"
        #                                               + "ex. /private/groups/kimlab/data/tempus/illumina/20241107/raw/metadata.csv"
        #                                               + "this file should not have a header, the first column should be sampleId, the second the sample type"
        # )

        # self.requiredArg.add_argument( '-g', '--genesOfInterest', required=True, nargs='+', metavar="",
        #                           action='store', 
        #                           help="list of genes of interest"
        # )

        self.requiredArg.add_argument( '-l', '--lfcThreshold', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=float,
                                            help="log fold change cut off"
        )  

        self.requiredArg.add_argument( '-n', '--number', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=int,
                                            help="the number of rows to be select"
        )  
        
        self.requiredArg.add_argument( '-p', '--padjThreshold', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=float,
                                            help="adjust p-value cut off"
        )   


###############################################################################
def main(inCommandLineArgsList=None):
    '''
    todo
    '''
    
    # we only configure logging in main module
    # loglevel = p.getProperty("LOG_LEVEL")
    loglevel = "INFO"
    # logFMT = p.getProperty("LOG_FMT")
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    logger = logging.getLogger(__file__)

    logger.warning("BEGIN")

    #
    # always log run time env to make debug easier
    #
    ORIG_PYTHONPATH = os.environ['PYTHONPATH']
    logger.warning(f'PYTHONPATH : {ORIG_PYTHONPATH}')
    logger.warning(f'FILE: {__file__}')
    logger.warning(f'PWD: {os.getcwd()}')

    cli = findBiomarkersCommandLine(
        version=__version__ , 
        author=__author__ ,
        date=__date__, 
        update=__updated__ )
    
    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    bioType       = cli.args.bioType
    countFilePath = cli.args.countFilePath
    padjThreshold = cli.args.padjThreshold
    lfcThreshold  = cli.args.lfcThreshold
    number        = cli.args.number  

    if not bioType is None and countFilePath is None:
        logger.error(f'--biotype requires --countFilePath')
        sys.exit(1)

    logger.warning(f'END')
    sys.exit(0)

################################################################################
if __name__ == '__main__':
    debugCommandLineArgsList=[
        #"--help",
        "--bioType", "protein_coding",
        "--padjThreshold", "0.001",
        "--lfcThreshold", "2.0",
        "--number" , "10",
    ]    
    main(debugCommandLineArgsList)

