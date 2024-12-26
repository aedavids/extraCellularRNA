#
# findBiomarkers.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref:

# findBiomarkers display the doc string 
'''
findBiomarkers: \n
select best biomarkers from a DESeq2 results file. The best biomarkers are selected based on the log fold change and adjusted p-value. The best biomarkers are written to biomarkerDESeq2Results.csv file.

'''

from analysis.bestSignatureGeneConfig import BestSignatureGeneConfig
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
                                                      required=False,
                                                      help="select biomarkers from biotype. Default is to select from, all genes"
        )        
        
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

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

        self.requiredArg.add_argument( '-d', '--deseq2ResultsFilePath', default=None, 
                                                        metavar="",
                                                        action='store', 
                                                        required=True, # only required if biotype is specified
                                                        help="path to a csv file with gene counts and bio type. "
                                                            + "ex. '/private/groups/kimlab/data/tempus/illumina/20241107/create/results/data/Control_vs_UD_deseq_results.csv'"
        )
 
        self.requiredArg.add_argument( '-o', '--outDir', default=None, 
                                                        metavar="",
                                                        action='store', 
                                                        required=True, # only required if biotype is specified
                                                        help="path to directory to write the biomarkerDESeq2Results.csv file to"
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

    bioType                 = cli.args.bioType
    padjThreshold           = cli.args.padjThreshold
    lfcThreshold            = cli.args.lfcThreshold
    number                  = cli.args.number  
    deseq2ResultsFilePath   = cli.args.deseq2ResultsFilePath
    outDir                  = cli.args.outDir


    # load the deseq2 results file
    deseqDF = pd.read_csv(deseq2ResultsFilePath)
    if 'gene_biotype' in deseqDF.columns:
        deseqDF['gene_biotype'] = deseqDF['gene_biotype'].astype('category')

    logger.info(f'deseqDF.head()\n{deseqDF.head()}')

    expectedBioTypes = ['Coding', 'DNA', 'LINE', 'LTR', 'Microsatellite', 'Other', 'SINE', 'lncRNA']
    if bioType is not expectedBioTypes:
        logger.warning(f'bioType: {bioType} is not in {expectedBioTypes}')

    #logger.info(f'deseqDF.loc[:, "gene_biotype"].cat.categories: {deseqDF.loc[:, "gene_biotype"].cat.categories}')

    if bioType is not None:
        deseqDF = deseqDF[ deseqDF['gene_biotype'] == bioType ]
        logger.info(f'deseqDF.head()\n{deseqDF.head()}')

    # set arguments we do not need to deprecated to make debugging easier
    # in the event we really need these arguments
    bsgc = BestSignatureGeneConfig(  
                            dataSetName="AEDWIP_deprecated", 
                            design="AEDWIP_deprecated", 
                            padjThreshold=padjThreshold, 
                            lfcThreshold=lfcThreshold,
                            n=number, 
                            localCacheRootPath="AEDWIP_deprecated", 
                            title="AEDWIP_deprecated"
        ) 

    biomarkerDF = bsgc.findGenes(deseqDF, "AEDWIP_deprecated")
    logger.info(f'best biomarkerDF\n{biomarkerDF}')

    os.makedirs(outDir, exist_ok=True)
    biomarkerPath = f'{outDir}/biomarkerDESeq2Results.csv'
    biomarkerDF.to_csv(biomarkerPath, index=False)

    logger.warning(f'END')
    sys.exit(0)

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

    main()

