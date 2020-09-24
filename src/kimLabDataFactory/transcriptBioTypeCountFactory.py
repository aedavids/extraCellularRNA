'''
Created on Sep 10, 2020

@author: andrewdavidson
'''
import pandas as pd

################################################################################
class TranscriptBioTypeCountFactory (object):
    '''
    ref:  juypterNotebooks/biotypeExploration.ipynb

    load count and bio type data, created by extraCellularRNA/R/notebooks/kras.ipsc_exo_vs_bulk.DESeq.Rmd

    example of input csvFile format
        "","HGNCT","HGNCG","BioType","TPM","NumReads"
        "1","DDX11L1-202","DDX11L1","lncRNA",0,0
        "2","DDX11L1-201","DDX11L1","transcribed_unprocessed_pseudogene",0,0
        "3","WASH7P-201","WASH7P","unprocessed_pseudogene",9.703864,249.755
        "4","MIR6859-1-201","MIR6859-1","miRNA",0,0
        
    public functions:
        
    '''

    ################################################################################
    def __init__(self, sampleType, fileList):
        '''
        arguments 
            sampleType
                a string
                example "bulk" or "exo"
            fileList:
                list of files to load class doc for example of expected file format
                for best results all files should be biologic replicants. For example
                
                ['../data/R/output/kras.ipsc.data.bulk.data.day.5.ctrl.1.biotype.csv',
                 '../data/R/output/kras.ipsc.data.bulk.data.day.5.ctrl.2.biotype.csv',
                 '../data/R/output/kras.ipsc.data.bulk.data.day.5.ctrl.3.biotype.csv']
        '''
        self.sampleType = sampleType
        self.fileList = fileList
        
    ################################################################################
    def loadFiles(self):
        """
        returns a pandas data frame. Use the 'replicant' column to identify row file source
        """
        dfList = self._loadFiles()
        DF = pd.concat(dfList, ignore_index=True)
        return DF
    
    #
    # private functions
    #
    
    ################################################################################    
    def _loadDataFrame(self, csvFile, replicant, numRows=None):
        """
        load count and bio type data, created by extraCellularRNA/R/notebooks/kras.ipsc_exo_vs_bulk.DESeq.Rmd
        
        arguments:
            csvFile:
            
            sampleType: string should be "bulk" or exo
            
            replicant: int. example 1, 2, or 3
            
        returns
            pandas data frame
            
        example of csvFile format
            "","HGNCT","HGNCG","BioType","TPM","NumReads"
            "1","DDX11L1-202","DDX11L1","lncRNA",0,0
            "2","DDX11L1-201","DDX11L1","transcribed_unprocessed_pseudogene",0,0
            "3","WASH7P-201","WASH7P","unprocessed_pseudogene",9.703864,249.755
            "4","MIR6859-1-201","MIR6859-1","miRNA",0,0
        """
        retDF = pd.read_csv(csvFile, nrows=numRows)
        # remove the "" column
        #print( [c for c in aedwipDF.columns] )
        colAxis = 1
        retDF = retDF.drop( ['Unnamed: 0'], axis=colAxis)
        retDF["sampleType"] = self.sampleType
        retDF["replicant"] = replicant
        
        # create factor variables
        retDF['BioType'] = retDF['BioType'].astype('category')
        retDF['sampleType'] = retDF['sampleType'].astype('category')
        #retDF['replicant'] = retDF['replicant'].astype('category')
                
        return retDF
    
    ################################################################################  
    def _loadFiles(self):  
        DFList = []
        for i in range(len(self.fileList)):
            replicantId = i + 1
            df = self._loadDataFrame(self.fileList[i], replicantId)
            DFList.append(df)
            
        return DFList