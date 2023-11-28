#
# cibersortMixtureFactory.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
# 10/4/23
#
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from kimLabUtils.baseCommandLine import BBaseCommandLine

###############################################################################
class CibersortSignatureMatrixFactoryCLI( BBaseCommandLine ):
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
        # self.parser.add_argument( '-p', '--prefix', default=None, metavar="", type=str, action="store", 
        #                           help="save output files with prefix")
        
        self.parser.add_argument( '-o', '--outDir', default=".", metavar="",
                                                      action='store', 
                                                      help="locaiton to write output files"
        )

    
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )
        
        self.requiredArg.add_argument( '-s', '--geneSignatureProfilesDataRootDir', required=True, default=None, metavar="",
                                                      action='store', 
                                                      help="the path the the geneSignatureProfiles files. "
                                                      + "There should be file for each category. "
                                                      + "The file contains the DESeq2 results for the genes "
                                                      + "of interest."
        )

        # self.requiredArg.add_argument( '-s', '--signatueGeneFilePath', required=True, default=None, metavar="",
        #                                               action='store', 
        #                                               help="path to a tsv file created by "
        #                                                     + " createCiberSortGeneSignatureMatrix.ipynb. "
        #                                                     + "Each row in this file coresponds to a 1vsAll DESeq "
        #                                                     + "result for a specific gene"
        #)

        self.requiredArg.add_argument( '-g', '--groupByGeneCountFilePath', required=True, default=None, metavar="",
                                                      action='store', 
                                                      help="path to a csv file with gene counts. "
                                                        + "ex. 'groupbyGeneTrainingSets/GTEx_TCGA_TrainGroupby.csv'"
                                                        + "first column is 'geneId' there is a column for each sample"
                                    )
        

        self.requiredArg.add_argument( '-c', '--colDataFilePath', required=True, default=None, metavar="",
                                                      action='store', 
                                                      help="path to a csv file containing sample meta data in DESeq format"
        )

        self.requiredArg.add_argument( '-f', '--scalingFactorsPath', required=True, default=None, metavar="",
                                                      action='store', 
                                                      help="path to a csv file of DESeq estimated scaling "
                                                      + "factors used to adjust each sample to account for "
                                                      + "libaray size and library composition"
        )

        self.requiredArg.add_argument( '-l', '--localCacheDir', required=True, default=None, metavar="",
                                                      action='store', 
                                                      help="when reading files will first check localCache. "
                                                      + "If not in cache will copy. Ensures file are localized."
        )

