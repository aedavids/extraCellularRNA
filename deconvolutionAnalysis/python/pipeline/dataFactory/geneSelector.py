#
# geneSelector.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
#

TODO delete this file

import logging
import pandas as pd
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

class GeneSelectorDEPRECATED( object ):
    '''
    Virtual Base class. Derived classes should implement findGenes(self, deseqDF : pd.DataFrame ) -> pd.DataFrame
    '''
    logger = logging.getLogger(__name__)

    ################################################################################    
    def __init__(self, ) :
        self.logger.debug(f'class does not have any data members to init')

    ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, 
                  signatureGeneConfig : SignatureGeneConfiguration ) -> pd.DataFrame:
        '''
        Find genes that that are statistically signifigant and up requlated
        
        arguments:
            deseqDF:
                results of DESeq2 as a pandas dataframe 
                
            signatureGeneConfig
                contains run parmeters            
                
        
        return:
            pandas dataframe
        '''
        self.logger.error("derived class did not implement findGenes()")
        raise NotImplementedError()
    