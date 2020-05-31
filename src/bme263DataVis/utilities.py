'''
Created on May 28, 2020

@author: andrewdavidson
'''

import matplotlib
import matplotlib.pyplot as plt
from os import path


###############################################################################
class MatPlotLibUtilities(object):
    '''
    classdocs
    '''

    ###############################################################################
    def __init__(self):
        '''
        Constructor
        '''
    ########################################################################
    def createPanel(self, fig,
                    panelWidthInInches, panelHeightInInches,
                    leftRelativeSize, bottomRelativeSize):
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
            [leftRelativeSize, bottomRelativeSize, relativeWidth, relativeHeight])
    
        return retPanel        
    
    ########################################################################
    def loadStyle(self):
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
        resourcesDir = path.join(path.dirname(__file__), 'styles')
        #print("AEDWIP resourcesDir:{}".format(resourcesDir))
        BME163MpltStylePath = path.join(resourcesDir, 'BME163.mplstyle')
        #print("AEDWIP BME163MpltStylePath:{}".format(BME163MpltStylePath))

        #plt.style.use('BME163.mpltstyle')        
        #plt.use.style(BME163MpltstylePath)
        plt.style.use(BME163MpltStylePath)
