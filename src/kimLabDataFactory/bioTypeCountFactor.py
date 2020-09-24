'''
Created on Sep 11, 2020

@author: andrewdavidson
'''

import pandas as pd

################################################################################
class BioTypeCountFactory(object):
    '''
    select pandas data frame rows in gene set
    
    ref:  
        juypterNotebooks/biotypeExploration.ipynb
        kimLabDataFactory.transcriptBioTypeCountFactory
        extraCellularRNA/R/notebooks/kras.ipsc_exo_vs_bulk.DESeq.Rmd

    public functions:
        __init__(self, csvFile)
        getGeneSet()
        selectRows(self, df)
    '''

    ################################################################################
    def __init__(self, csvFile):
        '''
        selects set of genes from csvFile
        
        arguments:
            csvFile
                example csv file:
                ../data/R/output/DESeq.ctrl.sampleType_ex0_vs_bulk.upRegulatedCounts.csv"
                
                (base) $ head DESeq.ctrl.sampleType_ex0_vs_bulk.upRegulatedCounts.csv
                "","bulk","bulk","bulk","ex0","ex0","ex0"
                "AC010970.1",4,5,2,65120,131914,61231
                "FP671120.6",10,9,7,258902,342833,316727
                "LINC01309",5,3,2,766,5584,1579
        '''
        self.geneSetSeries = None
        
        self._selectGeneSet(csvFile)
        
    ################################################################################
    def getGeneSet(self):
        '''
        returns pandas series
        '''
        return self.geneSetSeries

    ################################################################################
    def _selectGeneSet(self, csvFile):
        '''

        '''
        deseqDF = pd.read_csv(csvFile)

        # fix col names, change 'Unnamed: 0'
        origColNames = deseqDF.columns.values
        origColNames[0] = "HGNCG"
        deseqDF.columns = origColNames
        deseqHGNCGSeries = deseqDF.loc[:, "HGNCG"]
        self.geneSetSeries = deseqHGNCGSeries
        
    ################################################################################
    def selectRows(self, df):
        '''
        arguments:
            df: a pandas data frame loaded using the TranscriptBioTypeCountFactory class
            
        returns a pandas data frame
        
        '''
        selectRows = df.loc[:, "HGNCG"].isin(self.geneSetSeries)
        geneSetDF = df.loc[selectRows, :]
        return geneSetDF 
        
    
        
