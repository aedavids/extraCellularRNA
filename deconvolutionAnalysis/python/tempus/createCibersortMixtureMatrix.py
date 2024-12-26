#
# createCiberSortMixtureMatrix.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref:
# extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
# extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/cibersortMixtureMatrixFactory.py

# createCiberSortMixtureMatrices cli the doc string 
'''
create mixture.tsv file for cbersortx from count data created by 
"Create", formally known as 'Complete Seq'
    
ref:
    extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
    extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/cibersortMixtureMatrixFactory.py
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
__date__ = '2024-12-25'
__updated__ = '2024-12-25'

###############################################################################
class CreateMixtureMatrixCommandLine( BBaseCommandLine ):
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
        None

    
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

        self.requiredArg.add_argument( '-o', '--outDir', default=".", metavar="",
                                                      action='store', 
                                                      help="locaiton to write output file mixture.tsv file"
        )        

        self.requiredArg.add_argument( '-n', '--normalizedCountFilePath', required=True, default=None, metavar="",
                                                      action='store', 
                                                      help="path to a csv file with normalized gene counts. "
                                                        + "ex. '/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv'"
                                                        + "1st 3 column are gene_id, gene_name, gene_biotype there is a column for each sample"
        )

        self.requiredArg.add_argument( '-g', '--genesOfInterest', required=True, nargs='+', metavar="",
                                  action='store', 
                                  help="list of genes of interest"
        )

################################################################################
# def createMixtureMatrix(
#         countDF: pd.DataFrame,
#         genesOfInterest = list[str]):
#     '''
#     TODO
#     '''
#     pass

################################################################################
def main(inCommandLineArgsList=None):
    '''
    TODO
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

    #
    # parse the command line
    #
    cli = CreateMixtureMatrixCommandLine( 
        version=__version__,
        author=__author__, 
        date=__date__,
        update=__updated__ )

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    outDir = cli.args.outDir
    normalizedCountFilePath = cli.args.normalizedCountFilePath
    genesOfInterest = cli.args.genesOfInterest

    
    # make sure output director exist
    os.makedirs(outDir, exist_ok=True)

    #
    # load the normalized counts and get rid of any extra columns
    #
    normalizedCountDF = pd.read_csv(normalizedCountFilePath)

    normalizedCountDF.set_index("gene_id", inplace=True)
    
    if "gene_name" in normalizedCountDF.columns:
        normalizedCountDF.drop(columns=["gene_name"], inplace=True)

    if "gene_biotype" in normalizedCountDF.columns:
        normalizedCountDF.drop(columns=["gene_biotype"], inplace=True)

    logger.info(f'normalizedCountDF.shape : {normalizedCountDF.shape}')

    selectRows = normalizedCountDF.index.isin(genesOfInterest)
    retDF = normalizedCountDF.loc[selectRows, :]

    # cibersortx expects index name to be sampleTitle
    retDF.index.name = "sampleTitle"

    # sort the sample names to make analysis easier
    sortedCols = sorted(retDF.columns, reverse=True)
    retDF = retDF.loc[:, sortedCols]
    #
    # save output file
    #
    os.makedirs(outDir, exist_ok=True)
    outPath = f'{outDir}/mixtureMatrix.tsv'
    retDF.to_csv(outPath, index=True, sep='\t')

    print(f'saved mixture matrix to {outPath}')

    logger.warning("END")
    sys.exit(0)

################################################################################
if __name__ == '__main__':
    # debugCommandLineArgsList=[
    #     "--help",
    #     # "--outDir", "./tmp",
    #     # "--normalizedCountFilePath", "/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv",
    #     # "--genesOfInterest", "X7D_LINE", "Zaphod", 
    # ]    
    # main(debugCommandLineArgsList)

    main()
