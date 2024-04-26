#
# randomForestHyperparmeterSearch.py
# Andrew E. Davidson
# aedavids@ucsc.edu
# 04/07/2024
#

#
# calculateAUC()
# plotROCWrapper()
# plotROCWrapper()
#

import matplotlib.pyplot as plt
import numpy as np
import scikitplot as skplt
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

################################################################################
def calculateAUC ( 
        y : np.array,
        yProbability : np.array
    ) -> dict :
    '''
    ref: https://scikit-plot.readthedocs.io/en/stable/metrics.html#scikitplot.metrics.plot_roc

    arguments
        y : 
            ground truth labels
        
        yProbability : 
            predicted labels probablity  

    returns:
        dictionary :
            key = yProbability column idx
            value = area under ROC curve            
    '''

    retDict = dict()
    
    # print(f'yProbability.shape : {yProbability.shape}')
    numClasses = yProbability.shape[1]  # num columns
    for i in range(numClasses) :
        # The fpr, and tpr are only useful for plot the ROC curve
        # falsePositiveRate == specificity
        # truePositiveRate = sensitivity
        prob = yProbability[:,i]
        fpr, tpr, thresholds = roc_curve(y, prob)
        
        #area under the curve
        auc = roc_auc_score(y, prob)
        #print(f'\n************ calculateAUC() auc : {auc:.3f} ***********\n')
        
        retDict[i] = auc

    return retDict

################################################################################
def plotROC(
        panel : plt.axes, 
        y : np.array,
        yProbability : np.array,
        title : str,
        classesToPlot=None
    ) -> dict :
    '''
    plotsROC:

    ref: https://scikit-plot.readthedocs.io/en/stable/metrics.html#scikitplot.metrics.plot_roc

    TODO pass/inject function to ROCPlotFramework

    arguments:
        y : 
            ground truth labels
        
        yProbability : 
            predicted labels probablity

        classesToPlot :
            if binary classifier yProbability shape is (n,2)
            pass [1] to plot examples where y == 1

    returns:
        dictionary :
            key = yProbability column idx
            value = area under ROC curve
    '''
    rocPanel = skplt.metrics.plot_roc(y, yProbability,
                                    title=title,
                                    title_fontsize = "medium",
                                    text_fontsize = "small",
                                    ax = panel,
                                    plot_micro=False,
                                    plot_macro=False,
                                    classes_to_plot=classesToPlot
                                    )

    retDict = calculateAUC( y, yProbability )

    return retDict

################################################################################
def plotROCWrapper(XNP :np.array, 
                yNP :np.array, 
                model,
                title:str,
                 classesToPlot=None ) -> tuple[plt.figure, plt.axes, dict[int, float]]:
    '''
    arguments:
        XNP :
            samples to create prediction for
        y : 
            ground truth labels
    
        model:
            trained model 

    returns 
        fig:

        panel :
        
        dictionary :
            key = yProbability column idx
            value = area under ROC curve
 
    '''
    predictions  = model.predict(XNP)
    yProbability = model.predict_proba(XNP)

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(4, 4), ) #sharey=True
    
    #yProbability = yProbability[:,1]
    print(f'yNP.shape : {yNP.shape}')
    print(f'yProbability.shape : {yProbability.shape}' )
    
    classesToPlot = [1]
    aucDict = plotROC(
            axs, 
            yNP,
            yProbability,
            title,
            classesToPlot
    ) 

    return (fig, axs, aucDict)