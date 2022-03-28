'''
Created on Mar 27, 2022

@author: andrewdavidson
'''

__all__ = []
__version__ = 0.1
__date__ = '2022-03-26'
__updated__ = '2022-03-26'
__user_name__ = 'Andrew E. Davidson'

import pandas as pd
from plots.DESeqSelect import DESeqSelect
from plots.geneSignatureUpsetPlotCommandLine import GeneSignatureUpsetPlotCommandLine

########################################################################
def main( inComandLineArgsList=None ):
    '''
    process command line arguments load data and  call createPlot()
    '''
    cli = GeneSignatureUpsetPlotCommandLine( __user_name__, __version__, __date__, __updated__ )
    if inComandLineArgsList is None:
        cli.parse()
    else:
        cli.parse()( inComandLineArgsList )
        
    print(cli.args)
        
    numHeaderLines = cli.args.numHeaderLines
    files = cli.args.inputFiles
    
    for file in files:
        print("\nprocessing file: {}".format(file))
        tokens = file.split("/")
        #print(tokens)
        
        # last token: signatureGenesValidateThyroid.csv
        fileName = tokens[-1].split(".")[0]
        #print(fileName)
        tissueId = fileName[ len( "signatureGenesValidate" ):]
        print(tissueId)
        dataLoader = DESeqSelect( file )
        geneNamesNP, baseMeanNP, xlog2FoldChangeNP, yNeglog10pValueNP = dataLoader.readVolcanoPlotData(numHeaderLines)

########################################################################
if __name__ == '__main__':
    main()
