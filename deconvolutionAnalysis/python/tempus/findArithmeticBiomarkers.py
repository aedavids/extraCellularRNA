#
# findArithmeticBiomarkers.py
# Andrew E. Davidson
# aedavids@yucsc.edu
#

'''
findArithmeticBiomarkers.py\n

for tumor id: 
    finds the UD and Control samples from normalized counts
    saves (controlSeries - undilutedSeries).abs().sort_values(ascending=False).head(n=topN)
    to a file name arithmeticBiomarkers_{tumorId}.csv
'''

# from analysis.bestSignatureGeneConfig import BestSignatureGeneConfig
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
__date__ = '2025-01-23'
__updated__ = '2025-01-23'

###############################################################################
class findArithmeticBiomarkersCommandLine( BBaseCommandLine ):
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

        # self.parser.add_argument( '-b', '--bioType', default=None, metavar="",
        #                                               action='store', 
        #                                               required=False,
        #                                               help="select biomarkers from biotype. Default is to select from, all genes"
        # )        
        
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

        self.requiredArg.add_argument( '-c', '--normalizedCountsPath', default=None, 
                                                        metavar="",
                                                        action='store', 
                                                        required=True, 
                                                        help="path to a csv file with normalized gene counts and bio type. "
                                                            + "ex. '/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv'"
        )
 
       
        self.requiredArg.add_argument( '-o', '--outDir', default=None, 
                                                        metavar="",
                                                        action='store', 
                                                        required=True, # only required if biotype is specified
                                                        help="path to directory to write the arithmeticBiomarkers_{tumorId}.csv file to"
        )
        
        self.requiredArg.add_argument( '-n', '--topN', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=int,
                                            help="top n biomarkers to select"
        )  
 
        self.requiredArg.add_argument( '-t', '--tumorId', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=str,
                                            help="The tumor id. ex 'T1' or 'T2'"
        )  
          

 


###############################################################################
def main(inCommandLineArgsList=None):
    '''
    todo
    '''
    
    # we only configure logging in main module
    # loglevel = p.getProperty("LOG_LEVEL")
    loglevel = "INFO"
    # loglevel = "WARN"
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

    cli = findArithmeticBiomarkersCommandLine(
        version=__version__ , 
        author=__author__ ,
        date=__date__, 
        update=__updated__ )    

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    tumorId = cli.args.tumorId
    normalizedCountsPath = cli.args.normalizedCountsPath
    topN = cli.args.topN
    outDir = cli.args.outDir

    # start processing
    DF = pd.read_csv(normalizedCountsPath, index_col=0)
    logger.info(f'DF.shape: {DF.shape}')

    UDColName = DF.filter(like=tumorId + "_UD").columns[0]
    ControlColName = DF.filter(like=tumorId + "_Control").columns[0]
    logger.warning(f'UDColName: {UDColName}')
    logger.warning(f'ControlColName: {ControlColName}')

    controlSeries = DF.loc[:, ControlColName]
    undilutedSeries = DF.loc[:, UDColName]


    topBiomarkersSeries = (controlSeries - undilutedSeries).abs().sort_values(ascending=False).head(n=topN)
    
    topBiomarkersSeries.name = "arithmetic_diff"
    outPath = f'{outDir}/arithmeticBiomarkers_{tumorId}.csv'
    topBiomarkersSeries.to_csv(outPath, header=True)

    logger.warning(f'saved  {outPath} ')
    
    logger.warning("END")
    sys.exit(0)


################################################################################
if __name__ == '__main__':
    # debugCommandLineArgsList=[
    #     #"--help",
    #     "--normalizedCountsPath", "/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv",
    #     "--outDir", ".",
    #     "--topN", "3",
    #     "--tumorId", "T1"
    # ]    
    # main(debugCommandLineArgsList)

    main()
