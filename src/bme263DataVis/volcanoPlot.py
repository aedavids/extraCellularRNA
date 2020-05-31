#!/usr/local/bin/python3
# encoding: utf-8
'''
bme263DataVis.volcanoPlot.py -- shortdesc AEDWIP

It defines classes_and_methods

@author:     Andrew Davidson

@copyright:  2020 organization_name. All rights reserved.

@license:    license

@contact:    aedavids@ucsc.edu
@deffield    updated: Updated
'''
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from bme263DataVis.utilities import MatPlotLibUtilities
from kimLabDEQ.DESeqSelect import DESeqSelect

import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import os
import sys

__all__ = []
__version__ = 0.1
__date__ = '2020-05-26'
__updated__ = '2020-05-26'

###############################################################################
class CommandLine(object):
    '''
    Handle the command line, usage and help requests.
    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
    
        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''
    
        program_version = "v%s" % __version__
        program_build_date = str(__updated__)
        program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
        program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
        program_license = '''%s    
    
      Created by user_name on %s.
      Copyright 2020 organization_name. All rights reserved.
    
      Licensed under the Apache License 2.0
      http://www.apache.org/licenses/LICENSE-2.0
    
      Distributed on an "AS IS" basis without warranties
      or conditions of any kind, either express or implied.
    
    USAGE
    ''' % (program_shortdesc, str(__date__))    

        self.parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        self.parser.add_argument('-v', '--version', action='version', version=program_version_message)

        self.requiredArg = self.parser.add_argument_group('required arguments')
        
        # metavar
        # see https://stackoverflow.com/questions/26626799/pythons-argument-parser-printing-the-argument-name-in-upper-case
        self.requiredArg.add_argument( '-i', '--inputFile', required=True, default=None, metavar ="",
                                              action='store', help='input file name')
        self.requiredArg.add_argument('-o', '--outputFile', required=True, default=None, metavar ="",
                                             action='store', help='output file name')
    
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)
            
                
###############################################################################
class VolcanoPlot(object):
    '''
    classdocs
    '''

    ###############################################################################
    def __init__(self, panel):
        '''
        Constructor
        '''
        self.panel = panel
        
    ###############################################################################
    def plot(self, xlog2FoldChange, yNeglog10pValue):
        '''
        aedwip
        '''
        # we do not know what the true xlim should be
        # the values are not show in the template
        # this is a guess
#         self.panel.set_xlim(-12.5, 12.6)
#         self.panel.set_ylim(0, 60)

        # make sure all points are plotted
        xMin =  np.floor( np.min(xlog2FoldChange) )
        xMax = np.ceil( np.max(xlog2FoldChange) )
        self.panel.set_xlim(xMin, xMax)        
    
        # make sure all points are plotted
        yMin = np.floor( np.min(yNeglog10pValue) )
        yMax = np.ceil( np.max(yNeglog10pValue) )
        self.panel.set_ylim(yMin, yMax)        
            
        # https://matplotlib.org/tutorials/text/mathtext.html
        self.panel.set_xlabel(r'$log_2(fold\ change)$')
        self.panel.set_ylabel(r'$-log_{10}(adj\ p\ value)$')
    
        self.panel.plot(xlog2FoldChange,
                   yNeglog10pValue,
                   marker='o',
                   markerfacecolor='black',  # (56/255,66/255,156/255),
                   markeredgecolor='black',
                   markersize=1.5,  # diameter of mark # scatter is area
                   markeredgewidth=0,
                   linewidth=0)
    
#         self.panel.plot(extremeXlog2FoldChange,
#                    extremeYNeglog10pValue,
#                    marker='o',
#                    markerfacecolor='red',  # (56/255,66/255,156/255),
#                    markeredgecolor='red',
#                    markersize=1.5,  # diameter of mark # scatter is area
#                    markeredgewidth=0,
#                    linewidth=0)
#     
#         # add gene names to extreme  points
#         for i, geneName in enumerate(labels):
#             x = labelX[i]
#             y = labelY[i]
#             gn = geneName.strip() + " "
#             self.panel.text(x, y, gn, fontsize=6,
#                        horizontalalignment='right',
#                        verticalalignment='center')
#     
#         self.panel.tick_params(bottom=True, labelbottom=True,
#                           left=True, labelleft=True,
#                           right=False, labelright=False,
#                           top=False, labeltop=False)
    

########################################################################
def main(inComandLineArgsList=None):
    '''
    process command line arguments can call createPlot()  
    '''
    if inComandLineArgsList is None:
        cli = CommandLine()
    else:
        cli = CommandLine(inComandLineArgsList)
        

    mplu = MatPlotLibUtilities()
    mplu.loadStyle()

    inputFile = cli.args.inputFile
    dataLoader = DESeqSelect(inputFile)
    geneNamesNP, xlog2FoldChangeNP, yNeglog10pValueNP = dataLoader.readVolcanoPlotData()
    
    # standard paper size is 8.5 inches x 11 inches
    pageWidthInInches = 3
    pageHeightInInches = 3
    fig = plt.figure(figsize=(pageWidthInInches, pageHeightInInches))

    panelWidthInInches = 2
    panelHeightInInches = 2
    
    panel = mplu.createPanel(fig,
                        panelWidthInInches, panelHeightInInches,
                        leftRelativeSize=0.2, bottomRelativeSize=0.2)

    volcanoPlot = VolcanoPlot(panel)
    volcanoPlot.plot(xlog2FoldChangeNP, yNeglog10pValueNP)
    
    pre = "/public/groups/kimlab/"
    title = inputFile[len(pre):]
    panel.set_title(title, fontsize=8) #arial is not installed on courtyard, font is huge
    
    outputFile = cli.args.outputFile
    # png is an uncompressed bitmap format
    # output format is determined by output file name's suffix '.png'
    plt.savefig(outputFile)  # BME163 style sheet should set dpi=600

########################################################################
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    main()
