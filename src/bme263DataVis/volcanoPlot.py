#!/usr/local/bin/python3
# encoding: utf-8
'''
volcanoPlot.py -- shortdesc create a volcano plot

It defines classes_and_methods

@author:     Andrew Davidson

@copyright:  2020 organization_name. All rights reserved.

@license:    license

@contact:    aedavids@ucsc.edu
@deffield    updated: Updated
'''

from bme263DataVis.utilities import MatPlotLibUtilities
from bme263DataVis.volcanoPlotCommandLine import VolcanoPlotCommandLine
from kimLabDEQ.DESeqSelect import DESeqSelect

import matplotlib.pyplot as plt
import numpy as np

__all__ = []
__version__ = 0.1
__date__ = '2020-05-26'
__updated__ = '2020-05-26'
__user_name__ = 'Andrew E. Davidson'


###############################################################################
class VolcanoPlot( object ):
    '''
    classdocs
    '''

    ###############################################################################
    def __init__( self, mplu ):
        '''
        arguments:
            mplu: object of type MatPlotLibUtilities
            panel
        '''
        self.mplu = mplu
#         self.panel = panel

#     ###############################################################################
#     def createColorMapLedgend( self, colorBarPanel, minValueInDataCordinates,
#                          maxValueInDataCoordinates, numSteps, label ):
#         '''
#         AEDWIP
#         '''
#         # python issue, we can not pass a generator because
#         #  it is not indexable
#         yellow = tuple( ( c / 255 for c in [255, 255, 84] ) )
#         red = tuple( ( c / 255 for c in [255, 0, 0] ) )
#         green = tuple( ( c / 255 for c in [0, 255, 0] ) )
#         blueV = tuple( ( c / 255 for c in [  0, 9, 240] ) )
#         blue = tuple( ( c / 255 for c in [0, 0, 255] ) )
# 
#         colorMapTuple = self.mplu.createColorBar( colorBarPanel,
#                                                  minValueInDataCordinates,
#                                                  maxValueInDataCoordinates,
#                                                   yellow,
#                                                   red, numSteps, yLabel=label )
#         # RList, GList, BList = colorMapTuple
#         return colorMapTuple

    ###############################################################################
    def createColorMapLedgend( self, colorBarPanel,lowerColorTup, upperColorTup,
                                 minValueInDataCordinates, maxValueInDataCoordinates, numSteps, label ):
        '''
        AEDWIP
        '''
          
  
        colorMapTuple = self.mplu.createColorBar( colorBarPanel,
                                                 minValueInDataCordinates,
                                                 maxValueInDataCoordinates,
                                                  lowerColorTup,
                                                  upperColorTup, numSteps, yLabel=label )
        # RList, GList, BList = colorMapTuple
        return colorMapTuple
    
    ###############################################################################
    def plot( self, panel, xList, yList, colorByValues=None ):
        '''
        arguments:
            panel:
                type matplotlib.pyplot.axes

            x, y:
                array like object

            colorByValue:
                array type object of color tuples.  example color = (0.5, 0.3, 07)
                default = None
        '''
        # make sure all points are plotted
        xMin = np.floor( np.min( xList ) )
        xMax = np.ceil( np.max( xList ) )
        panel.set_xlim( xMin, xMax )

        # make sure all points are plotted
        yMin = np.floor( np.min( yList ) )
        yMax = np.ceil( np.max( yList ) )
        panel.set_ylim( yMin, yMax )

        # https://matplotlib.org/tutorials/text/mathtext.html
        panel.set_xlabel( r'$log_2(fold\ change)$' )
        panel.set_ylabel( r'$-log_{10}(adj\ p\ value)$' )

        if colorByValues is not None :
            # scatter markersize is area
            plotMarkersize = 0.6  # 0.75 #0.5 #1 #0.75 #0.5

            plotMarkerArea = np.pi * ( ( plotMarkersize / 2 ) ** 2 )
            scatterMarkersize = plotMarkerArea * 2  # strange did not look good with * 2
            panel.scatter( xList,
                          yList,
                          s=scatterMarkersize,
                          facecolor=colorByValues,
                          linewidth=0,
                          alpha=0.3 )

        else:
            # plot is faster than scatter how ever does not allow points to be
            # indvidually colored
            panel.plot( xList,
                       yList,
                       marker='o',
                       markerfacecolor='black',  # (56/255,66/255,156/255),
                       markeredgecolor='black',
                       markersize=1.5,  # diameter of mark
                       markeredgewidth=0,
                       linewidth=0,
                       alpha=0.3 )

#
#         self.panel.tick_params(bottom=True, labelbottom=True,
#                           left=True, labelleft=True,
#                           right=False, labelright=False,
#                           top=False, labeltop=False)


########################################################################
class _VolcanoPlotData :
    '''
    future proof private class use to manage data.

    avoid return list of parallel lists. Its hard to keep track of of kind of data is at
    each given position. using a class lets us name the various lists. Reduces refactor bugs
    '''

    def __init__( self ) :
        self.abundantX = []
        self.x = []

        self.abundantY = []
        self.y = []
        self.abundant = []

        self.abundantBaseMean = []
        self.baseMean = []

        # the xtick values for color gradient
        self.minAbundanceValue = 0
        self.maxAbundanceValue = 0


########################################################################
def loadData( inputFile ):
    '''
    returns a _VolcanoPlotData object
    '''
    ret = _VolcanoPlotData()

    dataLoader = DESeqSelect( inputFile )
    geneNamesNP, baseMeanNP, xlog2FoldChangeNP, yNeglog10pValueNP = dataLoader.readVolcanoPlotData()

    # find abundant values
    # create two data sets, values that are abundant and should be colored
    # and data that should be ploted as black points
    std = np.std( baseMeanNP )
    mean = np.mean( baseMeanNP )
    threshold = mean + 2 * std
    print( "AEDWIP mean:{} std:{} threshold:{}".format( mean, std, threshold ) )

    for i in range( len( baseMeanNP ) ):
        bm = baseMeanNP[i]
        x = xlog2FoldChangeNP[i]
        y = yNeglog10pValueNP[i]
        if bm >= threshold:
            ret.abundantX.append( x )
            ret.abundantY.append( y )
            ret.abundantBaseMean.append( bm )
        else:
            ret.x.append( x )
            ret.y.append( y )
            ret.baseMean.append( bm )

#     print( "AEDWIP mean:{} std:{} minCut:{}".format( mean, std, mean + 2 * std ) )

    # find the range of tick mark values for the color gradient panel
    ret.minAbundanceValue = int( np.floor( threshold ) )
    ret.maxAbundanceValue = int( np.ceil( np.max( baseMeanNP ) ) )

    return ret


########################################################################
def main( inComandLineArgsList=None ):
    '''
    process command line arguments load data and  call createPlot()
    '''
    cli = VolcanoPlotCommandLine( __user_name__, __version__, __date__, __updated__ )
    if inComandLineArgsList is None:
        cli.parse()
    else:
        cli.parse()( inComandLineArgsList )

    mplu = MatPlotLibUtilities()
    mplu.loadStyle()

    volcanoPlotData = loadData( cli.args.inputFile )

    # set up figure
    # standard paper size is 8.5 inches x 11 inches
    pageWidthInInches = 4
    pageHeightInInches = 3
    fig = plt.figure( figsize=( pageWidthInInches, pageHeightInInches ) )

    # create color map
    panelHeightInInches = 2
    bottomRelativeSize = 0.2
    leftRelativeSize = 0.9
    colorBarPanel = mplu.createPanel( fig,
                        3 / 16, panelHeightInInches,
                        leftRelativeSize, bottomRelativeSize )

    # set up the base mean color gradient legend
    label = "base mean"
    numSteps = 10  # 20  # 4  # 1000 # 4 #50 #20
    volcanoPlot = VolcanoPlot( mplu )
    
#         # python issue, we can not pass a generator because
#         #  it is not indexable
    yellow = tuple( ( c / 255 for c in [255, 255, 84] ) )
    red = tuple( ( c / 255 for c in [255, 0, 0] ) )    
    colorMapTuple = volcanoPlot.createColorMapLedgend( colorBarPanel,
                                                        yellow, red,
                                                        volcanoPlotData.minAbundanceValue,
                                                        volcanoPlotData.maxAbundanceValue,
                                                        numSteps, label )

    # create the color list data
    RList, GList, BList = colorMapTuple
    colorByValues = mplu.getColors( volcanoPlotData.abundantBaseMean,
                                   RList, GList, BList )

    black = ( 0, 0, 0 )
    allColors = [black] * len( volcanoPlotData.x ) + colorByValues

#     colorByValues = mplu.getColors( dataList, RList, GList, BList )

    # plot the main volcano plot
    panelWidthInInches = 2
    leftRelativeSize = 0.2
    volcanoPanel = mplu.createPanel( fig,
                        panelWidthInInches, panelHeightInInches,
                        leftRelativeSize, bottomRelativeSize )

    # put colored points at end of the list to avoid over plotting by black points
    allX = volcanoPlotData.x + volcanoPlotData.abundantX
    allY = volcanoPlotData.y + volcanoPlotData.abundantY
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
