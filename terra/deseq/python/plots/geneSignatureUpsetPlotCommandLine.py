
'''
plots.geneSignatureUpsetPlotCommandLine -- shortdesc

plots.geneSignatureUpsetPlotCommandLine is a description

It defines classes_and_methods

@author:     Andrew E. Davidson

@copyright:  2022 organization_name. All rights reserved.

@license:    license

@contact:    aedavids@ucsc.edu
@deffield    updated: Updated
'''


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from plots.volcanoPlotCommandLine import VolcanoPlotCommandLine

###############################################################################
class GeneSignatureUpsetPlotCommandLine( VolcanoPlotCommandLine ):
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
        super().__init__( author, version, date, update )
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
        # self.parser.add_argument( '-n', '--numHeaderLines', default=0, metavar="", type=int,
        #                                     action='store', help="number of header lines in input file")
        
        self.parser.add_argument( '-t', '--title', default=None, metavar="",
                                              action='store', help='plot title' )      
        
        self.parser.add_argument( '-w', '--width', default=10, type=int, action="store", help="image width in inches")
        
        self.parser.add_argument( '-z', '--height', default=4, type=int, action="store", metavar="HEIG", help="image height in inches")
        #  action="store",
        # default=3, type=int,
        # self.parser.add_argument( '-z', '--height', default=3, type=int, action="store", help="image height in inches")
        
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )
        
        self.requiredArg.add_argument('-d', '--dataSetsCSV', required=True, action='store',  metavar="",
                                      help="A csv file with 3 columns. the first column has setName,"
                                        + " The second contains the number of header lines to skip when reading."
                                        + " The third contains the path the the DESeqResults file")        
        
        # # good example for how to define vargs
        # # https://stackoverflow.com/a/60796254/4586180
        # self.requiredArg.add_argument('-i', '--inputFiles', required=True, action='store', nargs="+", 
        #                          help="one or more results files create by DESeq")
        
        # metavar
        # see https://stackoverflow.com/questions/26626799/pythons-argument-parser-printing-the-argument-name-in-upper-case
        # self.requiredArg.add_argument( '-i', '--inputFile', required=True, default=None, metavar="",
        #                                       action='store', help='input file name' )
        self.requiredArg.add_argument( '-o', '--outputFile', required=True, default=None, metavar="",
                                             action='store', help='output file name. suffix determines output format' )
        

        self.requiredArg.add_argument('-c', '--numThreads', required=True, default=None, type=int, metavar="", 
                                      action='store', help="currency level: size of thread pool")
        
        self.requiredArg.add_argument( '-i', '--intersectionOutputFile', required=True, default=None, metavar="",
                                       action='store', help="csv file containing intersection sets data")
