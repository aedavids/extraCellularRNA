#
# optimalSelectiveEnrichSignatureGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

from argparse import ArgumentParser
#https://stackoverflow.com/a/3853776/4586180
from argparse import RawDescriptionHelpFormatter
from kimLabUtils.baseCommandLine import BBaseCommandLine


###############################################################################
class OptimalSelectiveEnrichSignatureGeneConfigCLI( BBaseCommandLine ):
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

        # at least 1 or more var args. ie user defined args
        self.requiredArg.add_argument('-cl', '--classes', nargs="+",required=True, metavar="",
                                      action='store', 
                                      help="list of classes, types, or categories to enrich" \
                                            + "\n example: --classes UVM Whole_Blood"
                                      )
 
        self.requiredArg.add_argument( '-d', '--design', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="DESeq model design. Make sure the design uses docker friendly file path characters"
                                            +"\n Any alphanumeric characters from 0 to 9, A to Z, a to z, and the _ and - characters."
                                            + "\n example: 'tilda_gender_category"
        ) 

        self.requiredArg.add_argument( '-ds', '--dataSetName', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="data set name"
                                            + "\n example: GTEx_TCGA"
        )

        self.requiredArg.add_argument( '-hg', '--historicGeneSetPath', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="file containing a python set. use to ensure we do not add genes"
                                                + "\n that had been considered or removed by upstream stages"
                                            
        ) 

        self.requiredArg.add_argument( '-l', '--lfcThreshold', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=float,
                                            help="log fold change cut off"
        )  

        self.requiredArg.add_argument( '-m', '--maxNumberOfGenes', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=int,
                                            help="will try to add at most maxNumberOfGenes degree 1 genes"
        ) 
    
        self.requiredArg.add_argument( '-p', '--padjThreshold', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=float,
                                            help="adjust p-value cut off"
        )   

        self.requiredArg.add_argument( '-r', '--resultsDir', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="path to directory of deseq result files we will search for additional genes."
                                                + "\nThis should probably be the original results file. i.e. not the output of previous "
                                                + "\npipeline stages"
                                            
        ) 

        self.requiredArg.add_argument( '-s', '--startIdx', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=int,
                                            help="count start at zero. startIdx can be used to reduce the search space"
        ) 

        self.requiredArg.add_argument( '-t', '--title', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="upset plot title"
                                            + "\n example: a vs all"
        )
   
        self.requiredArg.add_argument( '-u', '--upStreamIntersectionDictionaryPath', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="output from upstream stage of pipeilne. "
                                                + "\nDictionary key is multi-index identifying sets that share elements. "
                                                + "\nvalue is the list of shared elements"
                                            
        ) 

        self.requiredArg.add_argument( '-w', '--windowLength', required=True, default=".", metavar="",   
                                            action='store', 
                                            type=int,
                                            help="defines the range of rows to search for candidate genes "
                                                + " \nThe range will be startIdx ... startIdx + windowLength"
        )  
          
