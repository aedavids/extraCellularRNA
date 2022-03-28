'''
Created on Mar 27, 2022

@author: andrewdavidson
'''

__all__ = []
__version__ = 0.1
__date__ = '2022-03-26'
__updated__ = '2022-03-26'
__user_name__ = 'Andrew E. Davidson'

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from plots.barChartCommandLine  import BarChartCommandLine
from plots.utilities import MatPlotLibUtilities

########################################################################
def main( inComandLineArgsList=None ):
    '''
    process command line arguments load data and  call createPlot()
    '''
    cli = BarChartCommandLine( __user_name__, __version__, __date__, __updated__ )
    if inComandLineArgsList is None:
        cli.parse()
    else:
        cli.parse()( inComandLineArgsList )
        
    print(cli.args)
    
    
    #
    # load data
    #
    dataDF = pd.read_csv( cli.args.inputFile )
    columnNames = dataDF.columns
    xLabel = columnNames[0]
    yLabel = columnNames[1]
    dataDF = dataDF.sort_values(yLabel)
    
    # xLabel will be something like lung, kidney, brain
    xTickLabels = dataDF.loc[:, xLabel].to_list()
    x = np.arange(len(xTickLabels))  + 1# the label locations

    
    # y values are the number of sample in each tissue type, i.e. int
    y = dataDF.loc[:,yLabel].to_list()
    print("y[0]:{}".format(y[0]))
        
    #
    # plot set up
    #
    mplu = MatPlotLibUtilities()
    mplu.loadStyle()
        
    figureWidthInInches = 8
    figureHeightInInches = 3
    fig = plt.figure(figsize=(figureWidthInInches,figureHeightInInches))   

    panelHeightInInches = figureHeightInInches #2
    bottomRelativeSize = 0.0 #0.2
    panelWidthInInches = figureWidthInInches #2
    leftRelativeSize = 0.0 #0.2
    barChartPanel = mplu.createPanel( fig,
                        panelWidthInInches, panelHeightInInches,
                        leftRelativeSize, bottomRelativeSize )
    
    #
    # config labels and axis
    #
    
    # 
    # xMin = np.floor(np.min( x ))
    # xMax = np.ceil( np.max(x) ) 
    # this mess up ax.bar() barChartPanel.set_xlim( xMin, xMax )    
    barChartPanel.set_xlabel( xLabel )
    barChartPanel.set_xticks(x) # , xLabels
    barChartPanel.set_xticklabels(xTickLabels, rotation = 90, fontsize=4)

    barChartPanel.set_ylabel( yLabel)
    # make sure all points are plotted
    # yMin = np.floor( np.min( y ) )
    # yMax = np.ceil( np.max( y ) )
    # this messes up ax.barChartPanel.set_ylim( yMin, yMax )
    
    #
    #
    # 
    width = 0.5 # default is width=0.8
    barChartPanel.bar(x, y, width=width) #, tick_label=xTickLabels
    # for i in range(len(x)):
    #     print("x:{} xTickLabels:{} y:{}".format(x[i], xTickLabels[i], y[i]))
            
    title = cli.args.title 
    if title:
        barChartPanel.set_title(title)
        

    outputFile = cli.args.outputFile
    fig.savefig(outputFile, dpi=300, bbox_inches='tight')
    print("saved plot: {}".format(outputFile))    
    
########################################################################
if __name__ == '__main__':
    main()
    