from argparse import ArgumentParser
#https://stackoverflow.com/a/3853776/4586180
from argparse import RawDescriptionHelpFormatter
from kimLabUtils.baseCommandLine import BBaseCommandLine


###############################################################################
class BestSignatureGeneConfigCommandLine( BBaseCommandLine ):
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

        # allow zero or more var args. ie user defined args
        #self.parser.add_argument('vargs', nargs="*", action='store')
    
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

        self.requiredArg.add_argument( '-c', '--localCacheRoot', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="local cache directory "
        ) 

        self.requiredArg.add_argument( '-d', '--design', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="DESeq model design. Make sure the design uses docker friendly file path characters"
                                            +"\n Any alphanumeric characters from 0 to 9, A to Z, a to z, and the _ and - characters."
                                            + "\n example: 'tilda_gender_category"
        ) 

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

        self.requiredArg.add_argument( '-s', '--dataSetName', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="data set name"
                                            + "\n example: GTEx_TCGA"
        )

        self.requiredArg.add_argument( '-t', '--title', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="upset plot title"
                                            + "\n example: a vs all"
        )
   
         
