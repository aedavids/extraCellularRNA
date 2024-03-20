#
# SignatureGeneConfiguration.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
#

import logging
import pandas as pd
from pipeline.dataFactory.utilities import urlify
import re
###############################################################################
class SignatureGeneConfiguration( object ):
    '''
    used to store all parameters needed to select sets of
    candidate biomarkers.

    should be treated as constant, immutable values

    public functions: see doc string for details
        findGenes(self, deseqDF : pd.DataFrame ) -> pd.DataFrame:

    '''

    logger = logging.getLogger(__name__)

    ################################################################################    
    def __init__(self, 
                 dataSetName : str, 
                 design : str, 
                 padjThreshold : float, 
                 lfcThreshold : float, 
                 n : int, 
                 localCacheRootPath : str, 
                 title : str
                ):
        '''
        arguments
            dataSetName: str
                
            Design:
                a string with the DESeq design. displayed on plots and encoded into data file names
                
            padjThreshold: float
                selects genes with padj values <= padjThreshold
                
            lfcThreshold : float
                log fold change threshold
            
            n: 
                type integer: 
                The number of rows to be select. 
                
            aedwip dataOutputBucketRoot
                example: gs://fc-e15b796f-1abe-4206-ab91-bd58374cc275/data/1vsAll/{up|down|best}
                location to store genes of interest
             ”   
            localCacheRootPath : str
                example: output
                
            title
                plot title       
                
        '''
        self.dataSetName = dataSetName
        self.design = design
        self.padjThreshold = padjThreshold
        self.lfcThreshold = lfcThreshold
        self.n = n
        self.localCacheRootPath = localCacheRootPath
        self.title = title
        
        # localCache = self.getLocalCachedDir()
        # ! mkdir -p $localCache

    ################################################################################
    def getfileNameBase(self) -> str:
        # docker can not mount file paths with special chars
        # https://www.ibm.com/docs/en/spectrum-symphony/7.3.1?topic=reference-docker-section
        #tmp = "{}-design:{}-padj:{}-lfc:{}-n:{}".format(
        tmp = "{} design {} padj {} lfc {} n {}".format(
                                self.dataSetName,
                                self.design,
                                self.padjThreshold,
                                self.lfcThreshold, 
                                self.n
                                )
        
        tmp = re.sub("~", " tilda ", tmp)
        tmp = re.sub("\\+", " ", tmp)

        
        tmp = urlify(tmp)
        return tmp

    # ################################################################################
    # def saveGenesOfInterestToBucketURL(self):
    #     return self.dataOutputBucketRoot + "/" + self.getfileNameBase() 

    ################################################################################
    def getLocalCachedDir(self) :
        return self.localCacheRootPath + "/" + self.getfileNameBase()                 

    ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame:
        '''
        Find genes that that are statistically and biologically signifigant 
        
        arguments:
            deseqDF:
                results of DESeq2 as a pandas dataframe          

           fileName:
                useful for logging, debugging, and addition down stream processing   
                
        
        return:
            pandas dataframe
        '''
        self.logger.error("derived class did not implement findGenes()")
        raise NotImplementedError()
   