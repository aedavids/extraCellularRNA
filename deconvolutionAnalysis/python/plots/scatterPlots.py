#
# scatterPlot.py
# Andrew E. Davidson
# aedavids@ucsc.edu
# 1/24/24
#

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

################################################################################
def createScatterPlot(
        x : list[float],
        y : list[float],
        title : str,
        xLabel : str,
        yLabel : str, 
        categories : list[str],
        #listOfCategoryTy : list[str],
        colors : list[str] = ['r',  'b'],
        figSize : tuple = (4,4),
        axisPaddingPrecentage : float = 0.02

    ) -> tuple[plt.figure, plt.axes]:
    '''
    '''

    fig = plt.figure(figsize=figSize)
    pannel = fig.add_subplot(1,1,1) 
    pannel.set_xlabel(xLabel, fontsize = 10)
    pannel.set_ylabel(yLabel, fontsize = 10)
    pannel.set_title(title, fontsize = 10)

    xMin = np.min(x)
    xMax = np.max(x)

    yMin = np.min(y)
    yMax = np.max(y)

    print(f'xmin: {xMin} xMax : {xMax}')
    print(f'yMin: {yMin} yMax : {yMax}')

    # move points off of edge of plot
    xBuffer = np.abs((xMax - xMin)) * axisPaddingPrecentage
    yBuffer = np.abs((yMax - yMin)) * axisPaddingPrecentage

    pannel.set_xlim(xMin - xBuffer, xMax + xBuffer)
    pannel.set_ylim(yMin - yBuffer, yMax + yBuffer)

    df = pd.DataFrame( {'x' : x, 
                        'y' : y,
                        'kind' : categories
                        })
    print(f'df.shape : {df.shape}')
    
    kinds = list( set(categories) )
    for kind, color in zip(kinds,colors):
        indicesToKeep = df['kind'] == kind
        pannel.scatter(
                    df.loc[indicesToKeep, 'x'],
                    df.loc[indicesToKeep, 'y'],
                    c = color,
                    #, s = 50
                    s=10, #5, #1,
                    linewidth=0 ,
                    alpha=0.5
                )
    pannel.legend(kinds)
    pannel.grid()

    return (fig, pannel)