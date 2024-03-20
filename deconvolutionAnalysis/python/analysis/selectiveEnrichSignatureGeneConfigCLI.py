#
# SelectiveEnrichSignatureGeneConfigCLI.py
# Andrew E. Davidson
# aedavids@uscs.edu
# 
from argparse import ArgumentParser
#https://stackoverflow.com/a/3853776/4586180
from argparse import RawDescriptionHelpFormatter
import logging
import os

from analysis.bestSignatureGeneConfigCLI import BestSignatureGeneConfigCommandLine

###############################################################################
class SelectiveEnrichSignatureGeneConfigCLI (BestSignatureGeneConfigCommandLine):
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

       # at least 1 or more var args. ie user defined args
        self.requiredArg.add_argument('-cl', '--classes', nargs="+",required=True, metavar="",
                                      action='store', 
                                      help="list of classes, types, or categories to enrich" \
                                            + "\n example: --classes UVM Whole_Blood"
                                      )
 
        self.requiredArg.add_argument( '-i', '--intersectionDict', required=True, default=".", metavar="",   
                                            action='store', 
                                            help="path to intersection dictionary file"
                                                + "\n see plots.test.testUpsetPlots testIntersection()"
                                                +"\n example: ./analysis/test/data/intersection.dict"
        ) 

        self.requiredArg.add_argument( '-a', '--numAdd', required=True, default="10", metavar="",   
                                            action='store', 
                                            type=int,
                                            help="integer. number of genes to add to each of the --classes"
        )          

    
