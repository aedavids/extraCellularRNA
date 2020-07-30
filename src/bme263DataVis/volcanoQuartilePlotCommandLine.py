'''
Created on Jun 17, 2020

@author: andrewdavidson
'''
from bme263DataVis.volcanoPlotCommandLine import VolcanoPlotCommandLine

###############################################################################
class VolcanoQuartilePlotCommandLine( VolcanoPlotCommandLine ):
    '''
    Handle the command line, usage and help requests.
    '''
    ###############################################################################
    def __init__( self, version, date, buildDate ):
        super().__init__(version, date, buildDate)

        
    ###############################################################################
    def _build( self ):
        super()._build()
        # added extra cli arguments here
