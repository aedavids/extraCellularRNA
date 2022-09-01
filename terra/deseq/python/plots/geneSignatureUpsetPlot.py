'''
Created on Mar 27, 2022

@author: andrewdavidson
'''
from pickle import TRUE

__all__ = []
__version__ = 0.1
__date__ = '2022-03-26'
__updated__ = '2022-03-26'
__user_name__ = 'Andrew E. Davidson'


#refactor see https://app.terra.bio/#workspaces/test-aedavids-proj/uber/notebooks/launch/signatureGenesUpsetPlots.ipynb?mode=edit

import logging
import os
import pandas as pd
from plots.DESeqSelect import DESeqSelect
from plots.geneSignatureUpsetPlotCommandLine import GeneSignatureUpsetPlotCommandLine
from matplotlib import pyplot as plt
import traceback
import upsetplot as up
from utils.upsetPlotData import UpSetPlotData

########################################################################
isloggingInitailized = False

def initializeLogging():
    '''
    insures logging will only be initialized once"
    '''
    print("CALLING initializeLogging()")
    global isloggingInitailized
    if not isloggingInitailized:
        print("first call initialize Logging")
        logging.basicConfig(level=logging.INFO, format="[%(levelname)s %(asctime)s %(filename)s:%(lineno)s - %(funcName) s()] %(message)s")
        isloggingInitailized = TRUE
        
initializeLogging()


########################################################################
class GeneSignatureUpsetPlot(object) :
    '''
    __init__(self, dataSetsDF, numThreads=1)
    
    plot(self, figureWidthInInches, figureHeightInInches)
        A easy of use function 
    
    getUpSetPlotData()
        returns object that can be use with the upsetplot python package's plot() function
        returns an object of type UpSetPlotData. 
        upData.plotData can be passed directly to the UpsetPlot plot function()
    
    getIntersectionDF()
        returns dataframe with information on sets and the genes
        columns: set name, gene name, DESeq Result values
        You may want to write this dataframe to a file

    getGeneSets()
        returns a dictionary [tissueId] = set( geneNames )
    
    '''
    logger = logging.getLogger(__name__)

    ################################################################################
    def __init__(self, dataSetsDF, numThreads=1):
        '''
        args:
            dataSetDF: type dataframe
                first column is the set name, second is file path the DESeq results
            
            numberThreads
                size of thread pool used to find set intersections
            
        '''
        initializeLogging()
        self.logger.info("BEGIN")
        
        self.dataSetsDF = dataSetsDF
        
        self.numThreads = numThreads
        
        # we want to be able to save the tissueIds, intersecting genes along with deseq values
        self.masterDataSet = None
        # the data we want to plot        
        self.geneSets = None
        self._createGeneSets()
        
        # transformed geneSet data into format upsetPlot requires
        self.upData = None
        self._createUpsetPlotData()
        
        # intersectionDF is dataframe with information on sets and the genes
        #  in their intersection
        self.intersectionDF = None
        self._createIntersectionDF()
        self.logger.info("END")
                    
    ################################################################################
    def getGeneSets(self):
        return self.geneSets
    
    ################################################################################
    def getIntersectionDF (self):
        '''
         returns dataframe with information on sets and the genes
         columns, set name, gene name, DESeq Result values
         '''
        return self.intersectionDF
    
    ################################################################################
    def getUpSetPlotData(self):
        '''
        returns object that can be use with the upsetplot python package's plot() function
        returns an object of type UpSetPlotData. 
        upData.plotData can be passed directly to the UpsetPlot plot function()
        '''
        return self.upData
    
    ################################################################################
    def plot(self, figureWidthInInches, figureHeightInInches):
        '''
        ref: https://upsetplot.readthedocs.io/en/stable/index.html
        
        returns tuple (figure, subplots)
            subplots:dict of matplotlib.axes.Axes
            Keys are ‘matrix’, ‘intersections’, ‘totals’, ‘shading’
        '''
        self.logger.info("BEGIN")
        fig = plt.figure(figsize=(figureWidthInInches,figureHeightInInches))
        subPlotDict = up.plot(self.upData.plotData, fig, element_size=None)
        self.logger.info("END")
        
        return (fig, subPlotDict)

        
    #
    # private functions
    #
    
    ################################################################################
    def _createGeneSets(self):
        self.masterDataSet = {} 
        self.geneSets = {} 
        numFiles = self.dataSetsDF.shape[0]
        
        for i in range(numFiles):
            tissueId = self.dataSetsDF.iloc[i,0]
            numHeaderLines = self.dataSetsDF.iloc[i,1]
            file = self.dataSetsDF.iloc[i,2]
            
            #tokens = file.split(".")
            isHack = False
            if "lfcShrink" in file:
                isHack = True
                
            self.logger.info("processing setId:{} numHeaderLines:{} file: {}".format(tissueId, numHeaderLines, file))
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
                
            self.geneSets[tissueId] = set( geneNamesNP )
            
            # hold on to deseq results Data so that we can analyze signature gene sets with overlapping genes
            self.masterDataSet[tissueId] = {
                "inputFile":file,
                # key is geneName
                "deseqResultSet": dataLoader.loadDESeqResultsAsStrings(numHeaderLines)
            }   
            
                 
    ################################################################################
    def _createIntersectionDF(self):
        '''
        TODO
        '''
        for setName, intersectionSet in self.upData.intersectionDict.items():
            tokens = setName.split(",")
            if len(tokens) > 1 :
                genesList = list(intersectionSet)
                # print("\n{}\n\t{}".format(setName, genesList))
                for geneName in genesList:
                    for tissueId in tokens:
                        deseqResult = self.masterDataSet[tissueId]['deseqResultSet']
                        stats = deseqResult[geneName]
                        # print("tissueId:{} gene:{} {}".format(tissueId, geneName, stats))
                        isHack = False
                        if "lfcShrink" in tissueId:
                            isHack = True
                            #print("tissueID:{}".format(tissueId))
                            
                        if isHack: 
                            # lfcShrink does not have a stat col
                            # print("ishack tissueID:{}".format(tissueId))
                            df = pd.DataFrame({
                                    "tissueId":[tissueId]
                                    ,"gene": [geneName]
                                     ,"baseMean": [stats[0]]
                                     ,"log2FoldChange": [stats[1]]
                                     ,"lfcSE": [stats[2]]
                                     ,"stat": "Na"
                                    ,"pvalue": [stats[3]]
                                    ,"padj": [stats[4]]
                                    })
                        
                        else :
                            # print("not isHack  tissueID:{}".format(tissueId))
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
                        if self.intersectionDF is not None :
                            self.intersectionDF = pd.concat( [self.intersectionDF, df] )
                        else:
                            self.intersectionDF = df        
        
    ################################################################################
    def _createUpsetPlotData(self):
        '''
        TODO AEDWIP
        '''
        try:
            self.upData = UpSetPlotData( self.geneSets, self.numThreads )
        except Exception as e:
            self.logger.error("UpSetPlotData failed")
            self.logger.error("exc:{}".format(e))
            self.logger.error(traceback.format_exc())
            
            raise e
        
########################################################################
def main( inComandLineArgsList=None ):
    '''
    process command line arguments load data and  call createPlot()
    '''
    logger = logging.getLogger(__name__)
   
    cli = GeneSignatureUpsetPlotCommandLine( __user_name__, __version__, __date__, __updated__ )
    if inComandLineArgsList is None:
        cli.parse()
    else:
        cli.parse()( inComandLineArgsList )
        
    logger.info("cli.args: {}".format(cli.args))
        
    
    # read data set CSV file
    # first column is the set name, second is file path the DESeq results
    logger.info("read data set CSV file: {}".format(cli.args.dataSetsCSV ))
    dataSetsDF = pd.read_csv( cli.args.dataSetsCSV )
    
    gsup = GeneSignatureUpsetPlot(dataSetsDF, cli.args.numThreads)

    logger.info("BEGIN plotting")
    figureWidthInInches = cli.args.width
    figureHeightInInches = cli.args.height
    # fig = plt.figure(figsize=(figureWidthInInches,figureHeightInInches))
    fig, subPlotDict = gsup.plot(figureWidthInInches, figureHeightInInches)
    
    logger.info("END plotting")
    
    
    title = cli.args.title
    if title:
        fig.suptitle( title, fontsize=8 )  # arial is not installed on courtyard, default font is huge

    outputFile = cli.args.outputFile
    fig.savefig(outputFile, dpi=300, bbox_inches='tight')
    print("saved plot: {}".format(outputFile))
    
    #
    # output sets that share genes and the genes in their intersection
    #
    logger.info("BEGIN creating intersection data set")

    logger.info("END creating intersection data set")
    intersectionOutputFile = cli.args.intersectionOutputFile
    intersectionDF = gsup.getIntersectionDF()
    intersectionDF.to_csv(intersectionOutputFile, index=False)
    print("\n*************** wrote file: {}".format(intersectionOutputFile)) 

########################################################################
if __name__ == '__main__':
    initializeLogging()
    main()
