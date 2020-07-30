#!/usr/local/bin/python3
# encoding: utf-8
'''
bme263DataVis.volcanoQuartilePlot.py -- shortdesc volcano plot where bottom and top quartiles are colored

It defines classes_and_methods

@author:     Andrew Davidson

@copyright:  2020 organization_name. All rights reserved.

@license:    license

@contact:    aedavids@ucsc.edu
@deffield    updated: Updated
'''

from bme263DataVis.utilities import MatPlotLibUtilities
from bme263DataVis.volcanoPlot import VolcanoPlot
from bme263DataVis.volcanoQuartilePlotCommandLine import VolcanoQuartilePlotCommandLine
from kimLabDEQ.DESeqSelect import DESeqSelect

import matplotlib.pyplot as plt
import numpy as np

__all__ = []
__version__ = 0.1
__date__ = '2020-06-18'
__updated__ = '2020-06-18'


########################################################################
class _VolcanoPlotData :
    '''
    future proof private class use to manage data.

    avoid return list of parallel lists. Its hard to keep track of of kind of data is at
    each given position. using a class lets us name the various lists. Reduces refactor bugs
    '''

    #######################################################################
    def __init__( self ) :
        self.lowerX = []
        self.upperX = []
        self.x = []

        self.lowerY = []
        self.upperY = []
        self.y = []

        self.lowerBaseMean = []
        self.upperBaseMean = []
        self.baseMean = []

        # the xtick values for color gradient
        self.minLowerValue = 0
        self.maxLowerValue = 0
        self.minUpperValue = 0
        self.maxUpperValue = 0

########################################################################
def createQtrColorMap(mplu, volcanoPlot, dataList, panel, minColor, maxColor, 
                      minValue, maxValue, label):
    # set up the base mean color gradient legend
    numSteps = 10  # 20  # 4  # 1000 # 4 #50 #20

    lowerColorMapTuple = volcanoPlot.createColorMapLedgend( panel,
                                                            minColor,
                                                            maxColor,
                                                            minValue,
                                                            maxValue,
                                                        numSteps, label )

    # create the color list data
    lowerRList, lowerGList, lowerBList = lowerColorMapTuple
    colorByValues = mplu.getColors( dataList, lowerRList, lowerGList, lowerBList )
    
    return colorByValues
    

########################################################################
def loadData( inputFile ):
    '''
    returns a _VolcanoPlotData object
    '''
    ret = _VolcanoPlotData()

    dataLoader = DESeqSelect( inputFile )
    geneNamesNP, baseMeanNP, xlog2FoldChangeNP, yNeglog10pValueNP = dataLoader.readVolcanoPlotData()

    lowerThreshold, upperThreshold = np.quantile( baseMeanNP, [0.25, 0.75 ] )

    for i in range( len( baseMeanNP ) ):
        bm = baseMeanNP[i]
        x = xlog2FoldChangeNP[i]
        y = yNeglog10pValueNP[i]
        if bm <= lowerThreshold:
            ret.lowerX.append( x )
            ret.lowerY.append( y )
            ret.lowerBaseMean.append( bm )
        elif bm >= upperThreshold:
            ret.upperX.append( x )
            ret.upperY.append( y )
            ret.upperBaseMean.append( bm )
        else:
            ret.x.append( x )
            ret.y.append( y )
            ret.baseMean.append( bm )

#     print( "AEDWIP mean:{} std:{} minCut:{}".format( mean, std, mean + 2 * std ) )

    # find the range of tick mark values for the color gradient panel
    ret.minLowerValue = 0
    ret.maxLowerValue = int( np.floor( lowerThreshold ) )
    ret.minUpperValue = int( np.floor( upperThreshold ) )
    ret.maxUpperValue = int( np.ceil( np.max( baseMeanNP ) ) )

    return ret


########################################################################
def main( inComandLineArgsList=None ):
    '''
    process command line arguments load data and  call createPlot()
    '''
    cli = VolcanoQuartilePlotCommandLine( __version__, __date__, __updated__ )
    if inComandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inComandLineArgsList )

    mplu = MatPlotLibUtilities()
    mplu.loadStyle()
    volcanoPlot = VolcanoPlot( mplu )

    volcanoPlotData = loadData( cli.args.inputFile )

    # set up figure
    # standard paper size is 8.5 inches x 11 inches
    pageWidthInInches = 4
    pageHeightInInches = 3
    fig = plt.figure( figsize=( pageWidthInInches, pageHeightInInches ) )

    #
    # lower quartile create color map
    #
    panelHeightInInches = 2
    bottomRelativeSize = 0.2
    leftRelativeSize = 0.09
    lowerQtrColorBarPanel = mplu.createPanel( fig,
                        3 / 16, panelHeightInInches,
                        leftRelativeSize, bottomRelativeSize )

    # python issue, we can not pass a generator because
    #  it is not indexable
    yellow = tuple( ( c / 255 for c in [255, 255, 84] ) )
    red = tuple( ( c / 255 for c in [255, 0, 0] ) )
    label = "lower quartile base mean"
    dataList = volcanoPlotData.lowerBaseMean
    lowerColorByValues = createQtrColorMap(mplu, volcanoPlot, dataList, lowerQtrColorBarPanel, red, yellow, 
                                           volcanoPlotData.minLowerValue, volcanoPlotData.maxLowerValue,
                                           label)
    #
    # upper quartile color map
    #
    leftRelativeSize = 0.91
    upperQtrColorBarPanel = mplu.createPanel( fig,
                        3 / 16, panelHeightInInches,
                        leftRelativeSize, bottomRelativeSize )
     
    label = "upper quartile base mean"
    dataList = volcanoPlotData.upperBaseMean
    green = tuple( ( c / 255 for c in [0, 255, 0] ) )
    blue = tuple( ( c / 255 for c in [0, 0, 255] ) )    
    upperColorByValues = createQtrColorMap(mplu, volcanoPlot, dataList, upperQtrColorBarPanel, blue, green, 
                                           volcanoPlotData.minUpperValue, volcanoPlotData.maxUpperValue,
                                           label)    


    
    #
    # plot the main volcano plot
    #
    black = ( 0, 0, 0 )
    allColors = [black] * len( volcanoPlotData.x ) + lowerColorByValues + upperColorByValues


    panelWidthInInches = 2
    leftRelativeSize = 0.25
    volcanoPanel = mplu.createPanel( fig,
                        panelWidthInInches, panelHeightInInches,
                        leftRelativeSize, bottomRelativeSize )

    # put colored points at end of the list to avoid over plotting by black points
    allX = volcanoPlotData.x + volcanoPlotData.lowerX + volcanoPlotData.upperX
    allY = volcanoPlotData.y + volcanoPlotData.lowerY + volcanoPlotData.upperY
    volcanoPlot.plot( volcanoPanel, allX, allY, allColors )

    title = cli.args.title
    if title:
        volcanoPanel.set_title( title, fontsize=8 )  # arial is not installed on courtyard, default font is huge

    # save image
    outputFile = cli.args.outputFile
    # png is an uncompressed bitmap format
    # output format is determined by output file name's suffix '.png'
    plt.savefig( outputFile, dpi=600 )  # BME163 style sheet should set dpi=600


########################################################################
if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    main()
