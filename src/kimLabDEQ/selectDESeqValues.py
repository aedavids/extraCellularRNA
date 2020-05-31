#!/usr/local/bin/python3
# encoding: utf-8
'''
kimLabDEQ.selectDESeqValues.py -- shortdesc AEDWIP

kimLabDEQ.select is a description AEDWIP

It defines classes_and_methods

@author:     Andrew Davidson

@copyright:  2020 organization_name. All rights reserved.

@license:    license

@contact:    aedavids@ucsc.edu
@deffield    updated: Updated
'''

import sys
import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import numpy as np
from kimLabDEQ.DESeqSelect import DESeqSelect

__all__ = []
__version__ = 0.1
__date__ = '2020-05-26'
__updated__ = '2020-05-26'

DEBUG = 0
TESTRUN = 0
PROFILE = 0

################################################################################
class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

################################################################################
def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2020 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        requiredArg = parser.add_argument_group('required arguments')

        # metavar
        # see https://stackoverflow.com/questions/26626799/pythons-argument-parser-printing-the-argument-name-in-upper-case
        requiredArg.add_argument("-a", "--adjPValue", dest="adjPValue", required=True ,
                            action="store", default=None, metavar ="",
                            help="select values with >= adjusted p-value")
        
        requiredArg.add_argument("-l", "--logFoldChange", required=True, 
                            dest="logFoldChange", action="store", default=None, metavar ="",
                            help="select values >= logFoldChange " )        
                
#         parser.add_argument("-r", "--recursive", dest="recurse", action="store_true", 
#                             help="recurse into subfolders [default: %(default)s]")
        
#         parser.add_argument("-v", "--verbose", dest="verbose", action="count",  default=0,
#                             help="set verbosity level [default: %(default)s]")
        
#         parser.add_argument("-i", "--include", dest="include", 
#                             help="only include paths matching this regex pattern. Note: exclude is given preference over include. [default: %(default)s]",
#                              metavar="RE" )
        
        '''parser.add_argument("-e", "--exclude", dest="exclude", 
                            help="exclude paths matching this regex pattern. [default: %(default)s]", 
                            metavar="RE" )'''
        
        parser.add_argument('-v', '--version', action='version', version=program_version_message)
        
        # metavar
        # see https://stackoverflow.com/questions/26626799/pythons-argument-parser-printing-the-argument-name-in-upper-case
        requiredArg.add_argument("-o", "--outputDir", dest="outputDir", required=True, 
                            action="store", default=None, metavar ="",
                            help="directory to create output files in.")
        
        requiredArg.add_argument(dest="paths", help="paths to one or more DESeq2 output files", 
                            metavar="path", nargs='+')

        # Process arguments
        args = parser.parse_args()

        paths = args.paths
        adjPValue = np.float(args.adjPValue)
        logFoldChange = np.float(args.logFoldChange)
         
        outputDir = args.outputDir
        for inputPath in paths:
            s = DESeqSelect(inputPath)
            s.saveRows(outputDir,  adjPValue, logFoldChange)
        return 0
    
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    
    except Exception as e:
        
        indent = len(program_name) * " "
        # AEDWIP TODO: FIXME: use logger
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help\n")
        raise(e)

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
        
    if TESTRUN:
        # https://docs.python.org/3/library/doctest.html
        import doctest
        doctest.testmod()
        
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'kimLabDEQ.select_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
        
    sys.exit(main())