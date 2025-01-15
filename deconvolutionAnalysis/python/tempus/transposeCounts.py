#
# transposeCounts.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref:

# trasposeCounts display the doc string 
'''
transposeCounts: \n
quick hack to reformat output of createTumorCountMatrix.sh into expected
DESeq2 count matrix format

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
__date__ = '2025-01-14'
__updated__ = '2025-01-14'

###############################################################################
class transposeCountsCommandLine( BBaseCommandLine ):
    '''
    Handle the command line, usage and help requests.
    '''
    ###############################################################################
    def __init__(self, version, author, date, update):
        '''
        Implement a parser to interpret the command line argv string using argparse.

        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''
        super().__init__(version, author, date, update )

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

        self.requiredArg.add_argument( '-c', '--createTumorCountMatrixFilePath', default=None, 
                                                        metavar="",
                                                        action='store', 
                                                        required=True, 
                                                        help="path to a csv file created by createTumorCountMatrix.sh. "
                                                            + "ex. '/private/groups/kimlab/data/tempus/illumina/20241107/create/results/data/Control_vs_UD_deseq_results.csv'"
        )
 
        self.requiredArg.add_argument( '-o', '--outFilePath', default=None, 
                                                        metavar="",
                                                        action='store', 
                                                        required=True, 
                                                        help="path to the transposed csv file to"
        )       
################################################################################
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

    cli = transposeCountsCommandLine(
        version=__version__ , 
        author=__author__ ,
        date=__date__, 
        update=__updated__ )

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')
    
    createTumorCountMatrixFilePath = cli.args.createTumorCountMatrixFilePath
    outFilePath = cli.args.outFilePath

    # head createTumorCountMatrix.sh.output/T1GroupByGenesCounts.csv |  cut -d , -f 1,2,3,4,76539,76540 > transposeTest.csv
    # tumorCountFilePath = "/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/tempus/bin/transposeTest.csv"
 
    tumorCountDF = pd.read_csv(createTumorCountMatrixFilePath,index_col="gene_id")
    logger.info(f'tumorCountDF.shape :{tumorCountDF.shape}')
    
    retDF = tumorCountDF.transpose()

    # geneId

    retDF.to_csv(outFilePath, index=True, index_label="geneId")
    logger.warning(f"saved transposed count matrix to : {outFilePath}")

    
    logger.warning("END")

################################################################################
if __name__ == '__main__':
    # debugCommandLineArgsList=[
    #     #"--help",
    #     "--bioType", "Coding",
    #     "--padjThreshold", "0.001",
    #     "--lfcThreshold", "2.0",
    #     "--number" , "10",
    #     "--deseq2ResultsFilePath", "/private/groups/kimlab/data/tempus/illumina/20241107/create/results/data/Control_vs_UD_deseq_results.csv"
    # ]    
    # main(debugCommandLineArgsList)

    debugCommandLineArgsList=[
        "--createTumorCountMatrixFilePath", "/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/tempus/bin/transposeTest.csv",
        "--outFilePath", "/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python/tempus/transposeTest.csv"
    ]
    main(debugCommandLineArgsList)

    #main()

