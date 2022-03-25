#!/usr/local/bin/python3
# encoding: utf-8

from   argparse import ArgumentParser
from   argparse import RawDescriptionHelpFormatter
import pandas as pd
import numpy as np

__all__ = []
__version__ = 0.1
__date__ = '2021-11-08'
__updated__ = '2021-11-08'

########################################################################
class CommandLine( object ):
    '''
    Handle the command line, usage and help requests.
    '''

    def __init__( self, inOpts=None ):
        '''
        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''

        program_version = "v%s" % __version__
        program_build_date = str( __updated__ )
        program_version_message = '%%(prog)s %s (%s)' % ( program_version, program_build_date )
        program_shortdesc = "creates a single count matrix file from a set of salmon quant files and calculates \n"\
        + "the DESeq2 equivalent scaling factors"
        # program_shortdesc = __import__( '__main__' ).__doc__.split( "\n" )[1]
        # print("WTF: {}".format(__import__( '__main__' ).__doc__.split( "\n" )))
        # print("AEDWIP program_shortdesc: " + program_shortdesc + "XXXX")

        program_license = '''%s

      Created by Andrew E. Davidson on %s.
      Copyright 2020 organization_name. All rights reserved.

      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0

      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.

    USAGE
    ''' % ( program_shortdesc, str( __date__ ) )

        self.parser = ArgumentParser( description=program_license, formatter_class=RawDescriptionHelpFormatter )
        self.parser.add_argument( '-v', '--version', action='version', version=program_version_message )

        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

        # metavar
        # see https://stackoverflow.com/questions/26626799/pythons-argument-parser-printing-the-argument-name-in-upper-case
        self.requiredArg.add_argument( '-t', '--transposeCSV', required=True, default=None, metavar="",
                                       action='store',
                                       help='file createSalmonNumReadsTransposeMatrix.sh' )

        self.requiredArg.add_argument( '-o', '--outputFile', required=True, default=None, metavar="",
                                       action='store', help='parent path directories must exist' )
        
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args( inOpts )

########################################################################
def main( inComandLineArgsList=None ):
    '''
    todo
    '''
    
    if inComandLineArgsList is None:
        # parse arguments from sys.argv
        cli = CommandLine()
    else:
        cli = CommandLine( inComandLineArgsList )

    print("INFO: arguments:\n {}\n".format(cli.args), flush=True)

    transposeCSV = cli.args.transposeCSV
    outputFile = cli.args.outputFile 
    
    #df = pd.read_csv( transposeCSV, index_col=False )
    # pandas will assign integer as the row index. given the names
    # are long creating index may be time consuming, let spark assign ints
    #df = pd.read_csv( transposeCSV )

    # load data frame in chunks
    # reduced memory overhead
    # make it possible to track progress
    # https://www.thiscodeworks.com/python-how-to-see-the-progress-bar-of-read_csv-stack-overflow-python/61b7fb4ebc0f1d0015ffd21e
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-chunking



    dtypeDict = {}
    dtypeDict[0] = str
    numRowsInQuantFile = 5387496 # these are columns in the input csv file
    for i in range(1, numRowsInQuantFile):
        dtypeDict[i] = np.float64

        
    ntrain = 10411 # number samples in training data set
    chunkSize = 5
    # pre allocate list length. append is slow
    chunkList = [None] * int(ntrain / chunkSize + 5) # add some fudge
    reportProgress = 10
    i = 0
    for chunk in pd.read_csv(transposeCSV, chunksize=chunkSize, dtype=dtypeDict):
        chunkList[i] = chunk
        if (i % reportProgress ) == 0:
            print("loaded chunk:{} chunkSize:{} \n".format(i, chunkSize), flush=True)
        i += 1


    print("INFO finished loading chunks i:{} chunkSize:{}".format(i, chunkSize), flush=True)
    df = pd.concat((f for f in chunkList), axis=0)
    nrows,ncols = df.shape    
    print("INFO loaded csv. numRows:{} numCols:{}\n".format(nrows, ncols), flush=True)
    chunkList = None
    print("INFO set chunkList=None to release memory\n", flush=True)

    
    tdf = df.transpose( copy=False )
    nrows,ncols = tdf.shape    
    print("INFO after transpose.  numRows:{} numCols:{}\n".format(nrows, ncols), flush=True)

    
    # #
    # # set column names = sample names and sort
    # # give the large number of sample keeping order
    # # sorted is the only way to ensure order match
    # # DESeq colData ordering
    # #
    
    # # tdf.columns = ['a', 'b', 'c', 'd']
    # # print("HACK tdf.columns: {}".format(tdf.columns) )

    # firstRow = tdf.iloc[0,:].values.tolist()
    # #print("HACK first row: {}".format(firstRow))
    # tdf.columns = firstRow
    # #print("HACK  set to first row tdf.columns: {}".format(tdf.columns) )

    # cols = tdf.columns
    # firstColName = cols[0]
    # sortedCols = [firstColName] + sorted( cols[1:] )
    # # print("HACK type(sortedCols):{}".format(sortedCols))

    # stdf = tdf.loc[:, sortedCols]
    # nrows,ncols = stdf.shape
    # print("INFO after sort.  numRows:{} numCols:{}\n".format(nrows, ncols), flush=True)    


    
    # after transpose the col names are 0,1,2, ...
    # the transcript names are the row indexes
    #stdf.to_csv(outputFile, header=False)

    tdf.to_csv(outputFile, header=False)    
    print("INFO output written to {}\n".format(outputFile), flush=True)
    print("INFO finished\n", flush=True)
        
    
########################################################################
if __name__ == "__main__":
    '''
    aedwip 0
    aedwip short discription
    '''
    
    # import sys;sys.argv = ['', 'Test.testName']
    main()
