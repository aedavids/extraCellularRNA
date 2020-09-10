'''
Created on Sep 9, 2020

@author: andrewdavidson aedavids@ucsc.edu
'''

import matplotlib.pyplot as plt
import numpy as np

################################################################################
class PlotProbMass(object):
    '''
    Plot the probability Mass Function
    
    public functions:
    __init__()
    
    '''

    ################################################################################
    def __init__(self, panel, DF, legendLabels=None, title=None):
        """
        creates a matplotlib line plot of the probability mass function
        
        TODO: replace with bar graph
        
        arguments:
            panel:
                type matplotlib axes. plot will be rendered in this panel
                
            DF: 
                pandas data frame
                Assumes row index is the states the random variable can assume and 
                the columns are separate distribution. E.G. control, treatement, ...
            
            legedLabels: 
                default None
                example [s[len("up"):] for s in upDataSets]
                
            title:
                type string
                default None
            
        """
        nRows, nCols = DF.shape
        x = np.arange( 1, nRows + 1, 1)

        dataSetLabels = DF.columns.to_numpy()
        panel.plot(x, DF, label= dataSetLabels)
        
        if (legendLabels != None) :
            panel.legend(legendLabels )
        else:
            panel.legend(dataSetLabels )

        ticks = [i for i in range(1, len(DF.index) + 1)]
        panel.set_xticks( ticks )
        panel.set_xticklabels(DF.index.to_numpy(), rotation = 45, ha="right")

        panel.set_ylabel("probability")

        panel.set_title( title )
        