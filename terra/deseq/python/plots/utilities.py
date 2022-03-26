'''
Created on May 28, 2020

@author: andrewdavidson
'''

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
from os import path


###############################################################################
class MatPlotLibUtilities( object ):
    '''
    classdocs
    '''

    ###############################################################################
    def __init__( self ):
        '''
        Constructor
        '''

    ########################################################################
    def createColorBar( self, panel, minY, maxY, firstColor, lastColor, numSteps, yLabel ):
        '''
        Creates color bar panel

        input
            panel: type axes

            minY, maxY: floats in data coordinates

            firstColor, lastColor.
                the color map gradient end values
                color example: yellow = tuple((c/255 for c in [255, 255,  84]))

        returns:
            color map tuples: (RedList, GreenList, BlueList)
        '''

        # the x axis tick values are not displayed
        # set them so that we know how to configure the rectangle width
        minX = 0
        maxX = 1
        panel.set_xlim( minX, maxX )

        panel.set_ylim( minY, maxY )

        tickList = np.linspace( minY, maxY, 10 + 1, dtype=np.int )
        panel.set_yticks( tickList )

        panel.set_ylabel( yLabel )

        # fil in color bar
        cMapTuple = self.createColorMap( firstColor, lastColor, numSteps )
        RList, GList, BList = cMapTuple

        # Rectangle parameters: left, bottom, width, height are in data coordinates
        # panel coordinates are relative to paper size
        width = maxX - minX
        height = ( maxY - minY ) / numSteps
        left = 0
        #print( "aedwip numSteps:{} w:{} h:{}".format( numSteps, width, height ) )
        for i in range( numSteps ):
            color = ( RList[i], GList[i], BList[i] )
            bottom = i * height + minY
            rect = mplpatches.Rectangle( ( left, bottom ), width, height,
                                facecolor=color,
                                edgecolor='black',
                                linewidth=0 )  
            panel.add_patch( rect )

        panel.tick_params( bottom=False, labelbottom=False,
                          left=True, labelleft=True,
                          right=False, labelright=False,
                          top=False, labeltop=False )

        return cMapTuple

    ########################################################################
    def createColorMap( self, startRGB, endRGB, numSteps ):
        '''
        input
            startRGB, endRBG: tuples
                example (104/255, 227/255, 251/255)

            numSteps: int

        return (R, G, B)
            where R, G, and B are tuples of length numSteps
        '''
        R = np.linspace( startRGB[0], endRGB[0], numSteps + 1 )
        G = np.linspace( startRGB[1], endRGB[1], numSteps + 1 )
        B = np.linspace( startRGB[2], endRGB[2], numSteps + 1 )

        return ( R, G, B )

    ########################################################################
    def createPanel( self, fig,
                    panelWidthInInches, panelHeightInInches,
                    leftRelativeSize, bottomRelativeSize ):
        '''
        returns a 'panel'. I.E. a graph component we can put stuff into
        do not use a plt.subplot. it is not flexible enough

        we can have multiple panels
        the values of left, bottom are relative to the size of the figure. they
        should be values between 0 and 1
        '''

        figWidth, figHeight = fig.get_size_inches()
        relativeWidth = panelWidthInInches / figWidth
        relativeHeight = panelHeightInInches / figHeight

        # left, bottom, width, and height are relative to size of figure,
        # they should be values in range [0,1]
        retPanel = plt.axes( 
            [leftRelativeSize, bottomRelativeSize, relativeWidth, relativeHeight] )

        return retPanel

    ########################################################################
    def createPanelSameSizeAsFig( self, fig):
        """
        create a panel the same size as the figure
        """
        figWidth, figHeight = fig.get_size_inches()
        panelWidthInInches = figWidth
        panelHeightInInches = figHeight
        leftRelativeSize = 0
        bottomRelativeSize = 1
        panel = self.createPanel(fig, panelWidthInInches, panelHeightInInches,
                leftRelativeSize, bottomRelativeSize)
        
        return panel
    
    ########################################################################
    def getColors( self, dataList, RList, GList, BList ):
        '''
        returns a list of colors.

        arguments:
            dataList:
                values in data coordinates

            RedList, GreenList, BlueList:
                color map gradient values
        returns
            colors: a list of colors
                example color = (0.5, 0.3, 07)
        '''
        colors = []
        numBins = len( RList )
        maxData = int( np.ceil( np.max( dataList ) ) )
        minData = int( np.floor( np.min( dataList ) ) )
        # print( "AEDWIP min:{} max:{}".format( minData, maxData ) )
        binSize = ( maxData - minData ) / numBins
        for data in dataList:
            adjVal = ( data - minData )
            colorIdx = int( adjVal / binSize )
            # print( "AEDWIP data:{:,} adjVal:{:,} binSize:{} colorIdx:{}"
            #       .format( data, adjVal, binSize, colorIdx ) )
            color = ( RList[colorIdx], GList[colorIdx], BList[colorIdx] )
            colors.append( color )

        # return np.array( colors )
        return colors

    ########################################################################
    def loadStyle( self ):
#         # https://stackoverflow.com/questions/6028000/how-to-read-a-static-file-from-inside-a-python-package
#         try:
#             import importlib.resources as pkg_resources
#         except ImportError:
#             # Try backported to PY<37 `importlib_resources`.
#             import importlib_resources as pkg_resources
#
#         from . import styles  # relative-import the *package* containing the templates
#
# #         template = pkg_resources.read_text(templates, 'temp_file')
# #         # or for a file-like stream:
# #         template = pkg_resources.open_text(templates, 'temp_file')
#

        # https://stackoverflow.com/questions/1011337/relative-file-paths-in-python-packages
        resourcesDir = path.join( path.dirname( __file__ ), 'styles' )
        # print("AEDWIP resourcesDir:{}".format(resourcesDir))
        BME163MpltStylePath = path.join( resourcesDir, 'BME163.mplstyle' )
        # print("AEDWIP BME163MpltStylePath:{}".format(BME163MpltStylePath))

        # plt.style.use('BME163.mpltstyle')
        # plt.use.style(BME163MpltstylePath)
        plt.style.use( BME163MpltStylePath )
        
    ########################################################################
    def getPanelDimensionInInches(self, fig, panel):
        """
        ref: https://stackoverflow.com/questions/19306510/determine-matplotlib-axis-size-in-pixels
        
        returns (width, height)
        """
        bbox = panel.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        
        return (width, height)
    
