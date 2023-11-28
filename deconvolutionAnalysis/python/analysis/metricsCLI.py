
# metricsCLI.py
# parse command line arguments needed to run fractions.CibersortResultsAsKWayClassifier
# 
# Andrew E. Davidson
# aedavids@ucsc.edu
# 11/13/23
#
from argparse import ArgumentParser
#https://stackoverflow.com/a/3853776/4586180
from argparse import RawDescriptionHelpFormatter
from kimLabUtils.baseCommandLine import BBaseCommandLine

###############################################################################
class MetricsCommandLine( BBaseCommandLine ):
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

    ###############################################################################
    def _build( self ):
        self.parser = ArgumentParser( description=self._getLicence(), formatter_class=RawDescriptionHelpFormatter )
    
        #
        # optional arguments
        #

        # self.parser.add_argument( '-p', '--prefix', default=None, metavar="", type=str, action="store", 
        #                           help="save output files with prefix")

    
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )
        
        self.requiredArg.add_argument( '-e', '--expectedFractionsPath', required=True, default=".", metavar="",
                                            action='store', 
                                            help="path to file containing expected fractions"
                                            + "\n example: analysis/test/data/expectedFractions.tsv"
        ) 

        self.requiredArg.add_argument( '-f', '--fractionsPath', required=True, default=".", metavar="",
                                            action='store', 
                                            help="path to CIBERSORTx fractions output"
                                            + "\n example: analysis/test/data/results.tsv"
        )     

        self.requiredArg.add_argument( '-o', '--outDir', required=True, default=".", metavar="",
                                            action='store', 
                                            help="locaiton to write output files"
        )
        