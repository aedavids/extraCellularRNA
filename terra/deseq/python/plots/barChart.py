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
    # print(dataDF)
    
    # xLabel will be something like lung, kidney, brain
    xLabels = dataDF.loc[:, xLabel].to_list()
    # print(xLabels)
    x = np.arange(len(xLabels))  # the label locations

    
    # y values are the number of sample in each tissue type, i.e. int
    y = dataDF.loc[:,yLabel].to_list()
    # print(y)
    
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
    xMin = 0
    xMax = np.ceil( len(x) )
    barChartPanel.set_xlim( xMin, xMax )    
    barChartPanel.set_xlabel( xLabel )
    barChartPanel.set_xticks(x) # , xLabels
    barChartPanel.set_xticklabels(xLabels, rotation = 90, fontsize=4)

    barChartPanel.set_ylabel( yLabel)
    # make sure all points are plotted
    yMin = np.floor( np.min( y ) )
    yMax = np.ceil( np.max( y ) )
    barChartPanel.set_ylim( yMin, yMax )
    
    #
    #
    # 
    numBars = len(y)
    barChartPanel.bar(x, y, tick_label=xLabels)
            
    title = cli.args.title 
    if title:
        barChartPanel.set_title(title)
        
    #barChartPanel.set_xticks(x, labels)

# x = np.arange(len(labels))  # the label locations
# width = 0.35  # the width of the bars
#
# fig, ax = plt.subplots()
# rects1 = ax.bar(x - width/2, men_means, width, label='Men')
# rects2 = ax.bar(x + width/2, women_means, width, label='Women')
#
# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Scores')
# ax.set_title('Scores by group and gender')
# ax.set_xticks(x, labels)
# ax.legend()
#
# ax.bar_label(rects1, padding=3)
# ax.bar_label(rects2, padding=3)

    outputFile = cli.args.outputFile
    fig.savefig(outputFile, dpi=300, bbox_inches='tight')
    print("saved plot: {}".format(outputFile))    
    
########################################################################
if __name__ == '__main__':
    main()
    