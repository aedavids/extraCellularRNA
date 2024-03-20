
#  exploreClassificationErrorsCLI.py
# parse command line arguments needed to run analysis.exploreClassificationErrors
# 
# Andrew E. Davidson
# aedavids@ucsc.edu
# 11/13/23
#
from argparse import ArgumentParser
#https://stackoverflow.com/a/3853776/4586180
from argparse import RawDescriptionHelpFormatter
from kimLabUtils.baseCommandLine import BBaseCommandLine

################################################################################
def _listOfStrings(args : str) -> list[str]:
    ret = args.split(",")
    return ret
  
###############################################################################
class ExploreClassificationErrorsCLI( BBaseCommandLine ):
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

        ref:
            https://docs.python.org/3/library/argparse.html 
            https://docs.python.org/3/library/argparse.html#sub-commands
        '''
        super().__init__(version, author, date, update )

    ###############################################################################
    def _build( self ):
        self.parser = ArgumentParser( description=self._getLicence(), formatter_class=RawDescriptionHelpFormatter )

        subparsers = self.parser.add_subparsers( dest="subcommand",
                                                title="subcommands",
                                                #description="valid subcommand",
                                                #help='fn: find false negatives, fp: find false positives, sg find shared genes'
                                                # dest="subcommand",
                                                required=True,
                                                #action='store'
                                                # dest="subcommand"
        )

        self._createParserForFalseNegative(subparsers)
        self._createParserForFalsePositives(subparsers)
        self._createParserForSharedGenes(subparsers)

     

    ################################################################################
    def _createParserForFalseNegative(self, subparsers ):
        parserFN = subparsers.add_parser('fn',
                                        # title='subcommands',
                                        # description='valid subcommands',
                                        help='find false negatives for class'
        )        

        self. _createSharedArguments(parserFN)
   

    ################################################################################
    def _createParserForFalsePositives(self, subparsers ):
        parserFP = subparsers.add_parser('fp',
                                        help='find false positives for class'
        )        

        self. _createSharedArguments(parserFP)

    ################################################################################
    def _createSharedArguments(self, parser):
        parser.add_argument( '-c', '--category',
                              required=True, 
                              default=None, 
                              metavar="",
                                action='store', 
                                help="The class, category, ... "
                                + "\n example: Whole_Blood"
        ) 

        parser.add_argument( '-e', '--expected',
                              required=True, 
                              default=None, 
                              metavar="",
                                action='store', 
                                help="path to extected fractions file"
                                + "\n example: analysis/test/data/expectedFractions.tsv"
        ) 

        parser.add_argument( '-o', '--outDir', required=True, default=".", metavar="",
                                            action='store', 
                                            help="locaiton to write output files"
        )

        parser.add_argument( '-r', '--results',
                                required=True, 
                                default=None, 
                                metavar="",
                                action='store', 
                                help="path to CIBERSORTx fractions result file "
                                + "\n example: analysis/test/data/results.tsv"
        )  
  
    ################################################################################
    def _createParserForSharedGenes(self, subparsers):
        sharedGenesParser = subparsers.add_parser('sg',
                                        help='find genes that are shared between a list of classes'
        )  

        # https://www.geeksforgeeks.org/how-to-pass-a-list-as-a-command-line-argument-with-argparse/
        
        sharedGenesParser.add_argument( '-c', '--category',
                                required=True, 
                                default=None, 
                                metavar="",
                                type=_listOfStrings,
                                action='store', 
                                help="a comma separated (no spaces) list of classes to search for. "
                                + "\n all interesections for which the list of classes is a subset will be returned"
                                + "\n example: -c Whole_Blood,UVM "
        )          

        sharedGenesParser.add_argument( '-i', '--intersectionDictionary',
                                required=True, 
                                default=None, 
                                metavar="",
                                action='store', 
                                help="path to file created by upsetplot that contains set interection information "
                                + "\n example: analysis/test/data/intersection.dict"
        )  

        sharedGenesParser.add_argument( '-o', '--outDir', required=True, default=".", metavar="",
                                            action='store', 
                                            help="locaiton to write output files"
        )
