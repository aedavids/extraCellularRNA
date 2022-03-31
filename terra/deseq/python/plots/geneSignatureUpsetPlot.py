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
        
    # numHeaderLines = cli.args.numHeaderLines
    
    # read data set CSV file
    # first column is the set name, second is file path the DESeq results
    dataSetsDF = pd.read_csv( cli.args.dataSetsCSV )
    print(dataSetsDF )
    
    
    masterDataSet = {} # we want to be able to save the tissueIds, intersecting genes along with deseq values
    geneSets = {} # the data we want to plot
    numFiles = dataSetsDF.shape[0]
    for i in range(numFiles):
        tissueId = dataSetsDF.iloc[i,0]
        numHeaderLines = dataSetsDF.iloc[i,1]
        file = dataSetsDF.iloc[i,2]
        
        #tokens = file.split(".")
        isHack = False
        if "lfcShrink" in file:
            isHack = True
            
        print("\nprocessing setId:{} numHeaderLines:{} file: {}".format(tissueId, numHeaderLines, file))
        # tokens = file.split("/")
        # #print(tokens)
        #
        # # last token: signatureGenesValidateThyroid.csv
        # fileName = tokens[-1].split(".")[0]
        # #print(fileName)
        # tissueId = fileName[ len( "signatureGenesValidate" ):]
        # print(tissueId)
        dataLoader = DESeqSelect( file )
        
        if isHack:
            geneNamesNP, baseMeanNP, xlog2FoldChangeNP, yNeglog10pValueNP = dataLoader.readVolcanoPlotData(numHeaderLines, hackPadjIndx=5)
        else:
            geneNamesNP, baseMeanNP, xlog2FoldChangeNP, yNeglog10pValueNP = dataLoader.readVolcanoPlotData(numHeaderLines)
            
        geneSets[tissueId] = set( geneNamesNP )
        
        # hold on to deseq results Data so that we can analyze signature gene sets with overlapping genes
        masterDataSet[tissueId] = {
            "inputFile":file,
            # key is geneName
            "deseqResultSet": dataLoader.loadDESeqResultsAsStrings(numHeaderLines)
            }
        
    upData = UpSetPlotData( geneSets )
    
    figureWidthInInches = 8
    figureHeightInInches = 3
    fig = plt.figure(figsize=(figureWidthInInches,figureHeightInInches))
    subPlotDict = up.plot(upData.plotData, fig)
    # [INFO] subPlotDict keys:dict_keys(['matrix', 'shading', 'totals', 'intersections'])

    # quick hack looked at plot , print names of intersecting genes
    print("\nThyroid intersection with Kidney")
    print( geneSets["Thyroid"].intersection(geneSets['Kidney_Cortex']) )
    
    print("\nThyroid intersection with Lung")
    print( geneSets["Thyroid"].intersection(geneSets['Lung']) )    
    
    title = cli.args.title
    if title:
        fig.suptitle( title, fontsize=8 )  # arial is not installed on courtyard, default font is huge

    outputFile = cli.args.outputFile
    fig.savefig(outputFile, dpi=300, bbox_inches='tight')
    print("saved plot: {}".format(outputFile))
    
    #
    # output sets that share genes and the genes in their intersection
    #
    print("\n\n\n $$$$$$$$$$$$$$$$ intersections:")
    intersectionDF = None
    for setName, intersectionSet in upData.intersectionDict.items():
        tokens = setName.split(",")
        if len(tokens) > 1 :
            genesList = list(intersectionSet)
            # print("\n{}\n\t{}".format(setName, genesList))
            for geneName in genesList:
                for tissueId in tokens:
                    deseqResult = masterDataSet[tissueId]['deseqResultSet']
                    stats = deseqResult[geneName]
                    # print("tissueId:{} gene:{} {}".format(tissueId, geneName, stats))
                    df = pd.DataFrame({
                            "tissueId":[tissueId]
                            ,"gene": [geneName]
                             ,"baseMean": [stats[0]]
                             ,"log2FoldChange": [stats[1]]
                             ,"lfcSE": [stats[2]]
                             ,"stat": [stats[3]]
                            ,"pvalue": [stats[4]]
                            ,"padj": [stats[5]]
                            })
                    if intersectionDF is not None :
                        intersectionDF = pd.concat( [intersectionDF, df] )
                    else:
                        intersectionDF = df
                        
    print( intersectionDF )
                                                                 

########################################################################
if __name__ == '__main__':
    main()
