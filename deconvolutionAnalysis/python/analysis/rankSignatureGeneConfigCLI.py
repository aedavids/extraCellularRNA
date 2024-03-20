from argparse import ArgumentParser
#https://stackoverflow.com/a/3853776/4586180
from argparse import RawDescriptionHelpFormatter
# from kimLabUtils.baseCommandLine import BBaseCommandLine
from analysis.bestRemoveHighDegreeSignatureGeneConfigCLI import BestSignatureGeneConfigCommandLine

###############################################################################
class RankSignatureGeneConfigCommandLine( BestSignatureGeneConfigCommandLine ):
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
        super()._build() 

        self.requiredArg.add_argument( '-r', '--deseqResultsDir', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="path to directory containing deseq results"
        )         

     
         
