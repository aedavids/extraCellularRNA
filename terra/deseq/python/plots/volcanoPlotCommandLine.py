'''
Created on Jun 17, 2020

@author: andrewdavidson
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter


###############################################################################
class VolcanoPlotCommandLine( object ):
    '''
    Handle the command line, usage and help requests.
    '''

    ###############################################################################
    def __init__( self, author, version, date, update ):
        '''
        Implement a parser to interpret the command line argv string using argparse.

        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''
        self.author = author
        self.program_version = version
        self.program_build_date = str( update )
        self.date = date

        self.program_version_message = '%%(prog)s %s (%s)' % ( self.program_version, self.program_build_date )
        self.program_shortdesc = __import__( '__main__' ).__doc__.split( "\n" )[1]

    ###############################################################################
    def _getLicence( self ):
        return '''%s

      Created by %s on %s.
      Copyright 2020 organization_name. All rights reserved.

      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0

      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.

    USAGE
    ''' % ( self.program_shortdesc, self.author, str(self.date ) )

    ###############################################################################
    def _build( self ):
        self.parser = ArgumentParser( description=self._getLicence(), formatter_class=RawDescriptionHelpFormatter )
                
        self.parser.add_argument( '-g', '--geneNamesFile',  default=None, metavar="",
                                              action='store', help='genes to color file. One gene name per line.' )
                
        self.parser.add_argument( '-n', '--numHeaderLines', default=0, metavar="", type=int,
                                            action='store', help="number of header lines in input file")
        
        self.parser.add_argument( '-t', '--title', default=None, metavar="",
                                              action='store', help='plot title' )
        
        self.parser.add_argument( '-v', '--version', action='version', version=self.program_version_message )
        
        # add a boolean argument
        # if present print gene name
        self.parser.add_argument( '-l', "--label", action="store_true", 
                                help="plot colored gene names") 

        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

        # metavar
        # see https://stackoverflow.com/questions/26626799/pythons-argument-parser-printing-the-argument-name-in-upper-case
        self.requiredArg.add_argument( '-i', '--inputFile', required=True, default=None, metavar="",
                                              action='store', help='input file name' )
        self.requiredArg.add_argument( '-o', '--outputFile', required=True, default=None, metavar="",
                                             action='store', help='output file name' )

    ###############################################################################
    def parse( self, inOpts=None ):
        self._build()
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args( inOpts )
