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
from plots.barChartCommandLine  import BarChartCommandLine
from matplotlib import pyplot as plt

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
    
    figureWidthInInches = 4
    figureHeightInInches = 3
    fig = plt.figure(figsize=(figureWidthInInches,figureHeightInInches))    
    
    #
    # load data
    #
    dataDF = pd.read_csv( cli.args.inputFile )
    columnNames = dataDF.columns
    # print(columnNames)
    
    xLabel = columnNames[0]
    # print(xLabel)
    x = dataDF.loc[:, xLabel]
    # print(x)
    
    yLabel = columnNames[1]
    # print(yLabel)
    y = dataDF.loc[:,yLabel]
    # print(y)
    
########################################################################
if __name__ == '__main__':
    main()
    