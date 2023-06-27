'''
Created on May 31, 2023

@author: andrewdavidson
aedavids@ucsc.edu

public functions:

def mapENSG_2_HUGO(
            txt2GeneFilePath : pl.Path,
            genesOfInterest=[] ,
        ) -> pd.DataFrame: 
        
def testMapENSG_2_HUGO():

def mapHUGO_2_ENSG(
            txt2GeneFilePath : pl.Path,
            genesOfInterest=[] ,
        ) -> pd.DataFrame:
'''

import pathlib as pl
import pandas as pd

################################################################################
def mapENSG_2_HUGO(
            txt2GeneFilePath : pl.Path,
            genesOfInterest=[] ,
        ) -> pd.DataFrame:
    '''
    
    arguments:
        txt2GeneFilePath : type pathLib Path or string
            2 columns, transcript name and gene id
            example 
            "/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv"
        
            ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|,DDX11L1
            
            hg38_rmsk_ALR/Alpha_range=chr22_KI270739v1_random:50071-73985_5'pad=0_3'pad=0_strand=-_repeatMasking=none,ALR/Alpha

        genesOfInterest
            type: list of strings
            default value = []
            
            
    example :
        genesOfInterest = ['ENSG00000223972.5', 'ENSG00000186081.12', 'ENSG00000205420.11'])
        
        returns data frame:
                                  ENSG     HUGO
            0        ENSG00000223972.5  DDX11L1
            136683  ENSG00000205420.11    KRT6A
            136688  ENSG00000186081.12     KRT5
    '''
    # countDF gene ids are in in ENSG format or loci . ex. ENSG00000227232.5 and (AAAAAC)n

    tx2GeneDF = pd.read_csv(txt2GeneFilePath, names=['uberId', 'geneId'])
    transcriptDF = tx2GeneDF

    # multiMapTranscriptDF = transcriptDF.loc[:, 'uberId'].str.split('|', expand=True)
    # multiMapTranscriptDF.head()
    # print("multi tail")
    # print(multiMapTranscriptDF.tail())

    # if you try and split a tx id that does not have "|" like the repeats
    # you get back columns filled with None
    mapDF = transcriptDF.loc[:, 'uberId'].str.split('|', expand=True).loc[:, [1, 5]]
    mapDF.columns = ['ENSG', 'HUGO']
    
    selectRows = mapDF.loc[:, "ENSG"].isin(genesOfInterest)
    HUGODF = mapDF.loc[selectRows, :]
    
    retDF = HUGODF.drop_duplicates(subset=['ENSG', 'HUGO'])
    return retDF

################################################################################
def testMapENSG_2_HUGO():
    tx2GenePath = "/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv"
    geneList = ['ENSG00000223972.5', 'ENSG00000186081.12', 'ENSG00000205420.11']
    hugoDF = mapENSG_2_HUGO(tx2GenePath, geneList)
    

    # TODO add assert
        # returns data frame:
        #                           ENSG     HUGO
        #     0        ENSG00000223972.5  DDX11L1
        #     136683  ENSG00000205420.11    KRT6A
        #     136688  ENSG00000186081.12     KRT5    
    print(hugoDF)

################################################################################
def mapHUGO_2_ENSG(
            txt2GeneFilePath : pl.Path,
            genesOfInterest=[] ,
        ) -> pd.DataFrame:
    '''
    
    TODO rework return. Make it same as MapENSG_2_HUGO. return DF that lets you map
    back and forth is more useful than a one way list. Prevents accidently re-odering bugs
    
    arguments:
        txt2GeneFilePath : type pathLib Path or string
            2 columns, transcript name and gene id
            example 
            "/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv"
        
            ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|,DDX11L1
            
            hg38_rmsk_ALR/Alpha_range=chr22_KI270739v1_random:50071-73985_5'pad=0_3'pad=0_strand=-_repeatMasking=none,ALR/Alpha

        genesOfInterest
            type: list of strings
            default value = []
        
    '''
    # countDF gene ids are in in ENSG format or loci . ex. ENSG00000227232.5 and (AAAAAC)n

    tx2GeneDF = pd.read_csv(txt2GeneFilePath, names=['uberId', 'geneId'])
    transcriptDF = tx2GeneDF
    # if len(genesOfInterest) > 0:
    #     transcriptDF = tx2GeneDF.loc[genesOfInterest, :]

    # multiMapTranscriptDF = transcriptDF.loc[:, 'uberId'].str.split('|', expand=True)
    # multiMapTranscriptDF.head()

    # if you try and split a tx id that does not have "|" like the repeats
    # you get back columns filled with None    
    mapDF = transcriptDF.loc[:, 'uberId'].str.split('|', expand=True).loc[:, [5, 1]]
    mapDF.columns = ['HUGO', 'ENSG']

    # ENSGDF.columns = ['ENSG']
    # # print(ENSGDF.shape)
    # # ENSGDF.head()
    #
    # ENSGList = list( ENSGDF.loc[:, "ENSG"].unique() )
    selectRows = mapDF.loc[:, "HUGO"].isin(genesOfInterest)
    ENSGDF = mapDF.loc[selectRows, :] 
    retDF = ENSGDF.drop_duplicates(subset=['HUGO', 'ENSG'])
    
    return retDF

################################################################################
if __name__ == '__main__':
    pass