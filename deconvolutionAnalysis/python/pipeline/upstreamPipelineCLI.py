#
# pipelineCLI.py
# parse command line arguments needed to run deconvolution hyper parameter tunning pipeline
# 
# Andrew E. Davidson
# aedavids@ucsc.edu
# 11/13/23
#
from argparse import ArgumentParser
#https://stackoverflow.com/a/3853776/4586180
from argparse import RawDescriptionHelpFormatter, REMAINDER
from kimLabUtils.baseCommandLine import BBaseCommandLine

###############################################################################
class UpstreamPipelineCommandLine( BBaseCommandLine ):
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

        self.parser.add_argument( '-u', '--useMedian', 
                                                      action='store_true', # if flag is present returns true else false
                                                      help="boolean, default = false. "
                                                      + "if true use mean, else use mean to calculate expected count values for each gene in a category"
        )

        # self.parser.add_argument( '-p', '--prefix', default=None, metavar="", type=str, action="store", 
        #                           help="save output files with prefix")

        # allow zero or more var args. ie user defined args
        # https://stackoverflow.com/a/25966881/4586180
        self.parser.add_argument('-v', '--vargs', nargs=REMAINDER, action='store')

        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )
        
        self.requiredArg.add_argument( '-c', '--colDataPath', required=True, default=".", metavar="",
                                            action='store', 
                                            help="path to deseq sample meta data file"
                                            + "\n example: pipeline/dataFactory/test/data/testIntegration/colData.csv "
        )     

        self.requiredArg.add_argument( '-s', '--countDataPath', required=True, default=".", metavar="",
                                            action='store', 
                                            help="path to count data file"
                                            + "\n example: pipeline/dataFactory/test/data/testIntegration/geneCounts.csv"
        ) 
     
        self.requiredArg.add_argument( '-d', '--deseqResultsDir', required=True, default=".", metavar="",
                                            action='store', 
                                            help="path to director containing 1vsAll results file."
                                            + "\n selects all files with names end with '.results' "
        )     
        self.requiredArg.add_argument( '-e', '--estimatedScalingFactors', required=True, default=".", metavar="",
                                            action='store', 
                                            help="path to director sample estimated scaling factors created by DESeq"
                                            + "\n example: pipeline/dataFactory/test/data/testSignatureGenes/1vsAll/estimatedSizeFactors.csv "
        )     

        self.requiredArg.add_argument( '-f', '--findModule', required=True, default=".", metavar="",
                                            action='store', 
                                            help="python module that implements createSignatureGeneConfig() "
                                            + "\nfunction. This function returns an object derived from "
                                            + "pipeline.dataFactory.signatureGeneConfig.  "
                                            + "\nYour class must implement findGenes(self, deseqDF : pd.DataFrame ) -> pd.DataFrame"
                                            + "\n"
                                            + "\nexample: --findModule pipeline.dataFactory.test.exampleCreateSignatureGeneConfig"
        )        

        self.requiredArg.add_argument( '-o', '--outDir', required=True, default=".", metavar="",
                                            action='store', 
                                            help="locaiton to write output files"
        )        

