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

import logging
import pathlib as pl
import pandas as pd

logger = logging.getLogger(__name__)

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
def mapHUGO_2_ENSG_DEPRECIATED(
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
    returns
        sample data frame
                    HUGO               ENSG                           bioType
            79436    PRELID1P1  ENSG00000217325.2              processed_transcript
            79437    PRELID1P1  ENSG00000217325.2  transcribed_processed_pseudogene
            157254     GOLGA8S  ENSG00000261739.2                   retained_intron
            157255     GOLGA8S  ENSG00000261739.2                    protein_coding
            179749     UBE2SP2  ENSG00000224126.2              processed_pseudogene
            195120  AC012615.3  ENSG00000267125.2                            lncRNA
            196526  AC010336.3  ENSG00000268120.1                            lncRNA
            227013     CCDC160  ENSG00000203952.9                    protein_coding]
    '''

    logger.info(f"BEGIN")
    logger.info(f'genesOfInterest :{genesOfInterest}')
    # countDF gene ids are in in ENSG format or loci . ex. ENSG00000227232.5 and (AAAAAC)n

    tx2GeneDF = pd.read_csv(txt2GeneFilePath, names=['uberId', 'geneId'])
    transcriptDF = tx2GeneDF
    # if len(genesOfInterest) > 0:
    #     transcriptDF = tx2GeneDF.loc[genesOfInterest, :]

    # multiMapTranscriptDF = transcriptDF.loc[:, 'uberId'].str.split('|', expand=True)
    # multiMapTranscriptDF.head()

    # if you try and split a tx id that does not have "|" like the repeats
    # you get back columns filled with None    
    geneCodeMapDF = transcriptDF.loc[:, 'uberId'].str.split('|', expand=True).loc[:, [5, 1, 7]]
    geneCodeMapDF.columns = ['HUGO', 'ENSG', 'bioType']
    geneCodeMapDF.dropna(inplace=True)
    geneCodeMapDF.drop_duplicates(inplace=True)

    # ENSGDF.columns = ['ENSG']
    # # print(ENSGDF.shape)
    # # ENSGDF.head()
    #
    # ENSGList = list( ENSGDF.loc[:, "ENSG"].unique() )

    if  len(genesOfInterest) > 0:
        selectRows = geneCodeMapDF.loc[:, "HUGO"].isin(genesOfInterest)
        ENSGDF = geneCodeMapDF.loc[selectRows, :] 
        retDF = ENSGDF.drop_duplicates(subset=['HUGO', 'ENSG'])

        logger.info(f'AEDWIP retDF:\n{retDF}')

        # log genes that did not map
        mappedGenesSeries = ENSGDF.loc[:, "HUGO"]
        logger.info(f'AEDWIP mappedGenesSeries :\n{mappedGenesSeries}')
        missingGenesRows = ~mappedGenesSeries.isin(genesOfInterest)
        logger.info(f'AEDWIP missingGenesRows :\n{missingGenesRows}')
        missingGenesSeries = mappedGenesSeries[missingGenesRows]
        logger.info(f'AEDWIP missingGenesSeries\n:{missingGenesSeries}')

        if missingGenesSeries.empty:
            logger.warning(f"!!!!!!! AEDWIP genes that did not map: {missingGenesSeries}")
    else:
        retDF = geneCodeMapDF
    
    logger.info(f"END")
    return retDF


################################################################################
def mapHUGO_2_ENSG(
            txt2GeneFilePath : pl.Path,
        ) -> pd.DataFrame:
    '''
    our original GTEx_TCGA biomarkers where quantified using HUGO name gencode.v35.ucsc.rmsk.tx.to.gene.csv.
    New verions of Complete-Seq using v39 and switched to ENSGO names. The function
    returns a dataframe that maps the HUGO names to something that can be mapped to v39

    The function is poorly name. The mapping files have two section. a top section of
    genecode anotations and a bottom section that comes the ucsc repeat mask. The formats are different
    The repeat mask gene name are not hugo or ENSGO names.

    Notice the repeat mask rows have 2 columns. The value in the second column is returns
    under the "HUGO" and "ENGS" columns even though it is not an idea from either of these 
    references. This makes find the mapping between gencode.v35.ucsc.rmsk.tx.to.gene.csv and v39 easier 
    
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

    TODO: check that we do not loose rows. Use grep to find the a couple of rows above
    and below the first repeat

    returns
        sample data frame
                        HUGO               ENSG                           bioType
        79436      PRELID1P1  ENSG00000217325.2              processed_transcript
        79437      PRELID1P1  ENSG00000217325.2  transcribed_processed_pseudogene
        157254       GOLGA8S  ENSG00000261739.2                   retained_intron
        157255       GOLGA8S  ENSG00000261739.2                    protein_coding
        179749       UBE2SP2  ENSG00000224126.2              processed_pseudogene
        195120    AC012615.3  ENSG00000267125.2                            lncRNA
        196526    AC010336.3  ENSG00000268120.1                            lncRNA
        227013       CCDC160  ENSG00000203952.9                    protein_coding
        230512         (TA)n              (TA)n                        repeatMask
        234711        LTR106             LTR106                        repeatMask
        238552         MER5C              MER5C                        repeatMask
        245986  HERVFH19-int       HERVFH19-int                        repeatMask
    '''

    logger.info(f"BEGIN")

    tx2GeneDF = pd.read_csv(txt2GeneFilePath, names=['uberId', 'geneId'])

    # if you try and split a tx id that does not have "|" like the repeats
    # you get back columns filled with None    
    # we can use this to separate the genecode and the repeat mask porition of the data frame
    geneCodeMapDF = tx2GeneDF.loc[:, 'uberId'].str.split('|', expand=True).loc[:, [5, 1, 7]]
    geneCodeMapDF.columns = ['HUGO', 'ENSG', 'bioType']

    # if HUGO, ENSG, and bioType == None, the row is part of the repeat mask
    selectRepeatMapRows = (geneCodeMapDF.loc[:,"HUGO"].isna()) & \
                            (geneCodeMapDF.loc[:,"ENSG"].isna()) & \
                            (geneCodeMapDF.loc[:,"bioType"].isna() )

    reapeatMaskDF = tx2GeneDF.loc[selectRepeatMapRows, :].copy()

    # clean up geneCodeMapDF. drop the repeat mask porition
    geneCodeMapDF.dropna(inplace=True)
    geneCodeMapDF.drop_duplicates(inplace=True)

    logger.info(f'tx2GeneDF.shape : {tx2GeneDF.shape}')
    logger.info(f'geneCodeMapDF.shape : {geneCodeMapDF.shape}')
    logger.info(f'repeatMaskDF.shape : {reapeatMaskDF.shape}')
   
    #
    # the genecode rows have three columns
    # the repeat mask only has 2
    # make the repeat mask portion look like the genecode dataframe
    # by renaming the geneId column to "HUGO", create a copy of the HUGO
    # column named ENSG. Create a column of "repeatMask" named 'bioType'
    #
    reapeatMaskDF["HUGO"] = reapeatMaskDF['geneId']
    nRows = reapeatMaskDF.shape[0]
    bioTypes = ["repeatMask"]*nRows
    reapeatMaskDF.loc[:, "bioType"] = bioTypes

    reapeatMaskDF = reapeatMaskDF.rename( columns={'geneId' : 'ENSG'} )
    reapeatMaskDF = reapeatMaskDF.drop( columns=['uberId'] )

    # reorder the columns to match geneCode dataframe
    reapeatMaskDF = reapeatMaskDF.loc[:, ['HUGO', 'ENSG', 'bioType']]

    retDF = pd.concat( [geneCodeMapDF, reapeatMaskDF] )

    logger.info(f"END")
    return retDF

################################################################################
if __name__ == '__main__':
    pass