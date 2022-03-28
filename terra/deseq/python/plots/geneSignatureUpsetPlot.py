'''
Created on Mar 27, 2022

@author: andrewdavidson
'''

__all__ = []
__version__ = 0.1
__date__ = '2022-03-26'
__updated__ = '2022-03-26'
__user_name__ = 'Andrew E. Davidson'

from plots.DESeqSelect import DESeqSelect
from plots.geneSignatureUpsetPlotCommandLine import GeneSignatureUpsetPlotCommandLine
from matplotlib import pyplot as plt
import upsetplot as up
from utils.upsetPlotData import UpSetPlotData

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
    
    geneSets = {}
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
        geneSets[tissueId] = set( geneNamesNP )
        
    upData = UpSetPlotData( geneSets )
    
    figureWidthInInches = 8
    figureHeightInInches = 3
    fig = plt.figure(figsize=(figureWidthInInches,figureHeightInInches))
    subPlotDict = up.plot(upData.plotData, fig)
    # [INFO] subPlotDict keys:dict_keys(['matrix', 'shading', 'totals', 'intersections'])

    # quick hack looked at plot , print names of intersecting genes
    print("\nThyroid intersection with Kidney_Cortex")
    print( geneSets["Thyroid"].intersection(geneSets['Kidney_Cortex']) )
    
    print("\nThyroid intersection with Lung")
    print( geneSets["Thyroid"].intersection(geneSets['Lung']) )    
    
    outputFile = cli.args.outputFile
    fig.savefig(outputFile, dpi=300, bbox_inches='tight')
    print("saved plot: {}".format(outputFile))

########################################################################
if __name__ == '__main__':
    main()
