#
# plasmaGAN.py
# Andrew E. Davidson
# aedavids@ucsc.edu
#
'''
public functions 
    loadCountData()
   
'''

import logging
import os.path
import pandas as pd

from analysis.utilities import loadList
from analysis.utilities import saveList

from intraExtraRNA.elifeUtilities import loadCounts
from intraExtraRNA.elifeUtilities import loadMetaData
from intraExtraRNA.elifeUtilities import selectFeatures

logger = logging.getLogger(__name__)

################################################################################
def loadCountData(
    localCacheDir : str,
    HUGO_featuresNames : list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str], pd.DataFrame]:
    '''
    given a list of HUGO gene biomarkers loads corresponding elife normalized
    count data.

    arguments:
        localCacheDir : str
            if not none
                if not in local cache 
                    load from /private and store in cache

        HUGO_featuresNames : list[str]

    returns :
        (XDF, metaDF, elifeLungGenes, missingElifeGenes, mapDF)

        XDF : index = sample id, columns = [gene_loci_name], values = counts

        metaDF : columns = [sample_id, diagnosis]

        elifeLungGenes : list of ensembl gene names

        missingElifeGenes

        mapDF : columns = [HUGO_v35	ENSG_v35	ENSG_v39]
    '''
    logger.info( "BEGIN" )

    XDFPath               = os.path.join( localCacheDir, "XDF.csv" )
    metaDFPath            = os.path.join( localCacheDir, "metaDF.csv" )
    elifeLungGenesPath    = os.path.join( localCacheDir, "elifeLungGenes.txt" )
    missingElifeGenesPath = os.path.join( localCacheDir, "missingElifeGenes.txt" )
    mapDFPath             = os.path.join( localCacheDir, "mapDF.csv" )

    if localCacheDir == None :
        logger.info( f'localCacheDir undefined loading from private')
        retTuple = loadCountDataImpl( HUGO_featuresNames )

        # if we do not set index. name pca.fit_transform fails
        XDF, metaDF, elifeLungGenes, missingElifeGenes, mapDF
        XDF.index.name = "sampleId"


    else :
        if os.path.isfile(XDFPath):
            logger.info( f'loading from {localCacheDir}' )
            # read everything from local cache
            XDF = pd.read_csv( XDFPath, index_col=0 )
            metaDF = pd.read_csv( metaDFPath )
            mapDF = pd.read_csv( mapDFPath )
            
            elifeLungGenes = loadList( elifeLungGenesPath )
            missingElifeGenes = loadList( missingElifeGenesPath )

            retTuple = (XDF, metaDF, elifeLungGenes, missingElifeGenes, mapDF)

        else :
            # load from private and store in cache
            logger.info( f'localCacheDir is empty, load from private')
            retTuple = loadCountDataImpl( HUGO_featuresNames )
            XDF, metaDF, elifeLungGenes, missingElifeGenes, mapDF = retTuple

            # save
            os.makedirs( localCacheDir, exist_ok=True )

            XDF.index.name = "sampleId"
            XDF.to_csv( XDFPath)
            metaDF.to_csv( metaDFPath )

            saveList( elifeLungGenesPath, elifeLungGenes, isSingleItemLine=False )
            saveList( missingElifeGenesPath, missingElifeGenes, isSingleItemLine=False )

            mapDF.to_csv( mapDFPath )


    logger.info( "END" )
    return retTuple

################################################################################
def loadCountDataImpl(
    HUGO_featuresNames : list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str], pd.DataFrame]:
    '''
    given a list of HUGO gene biomarkers loads corresponding elife normalized
    count data.

    arguments:
        HUGO_featuresNames : list[str]

    returns :
        (XDF, metaDF, elifeLungGenes, missingElifeGenes, mapDF)

        XDF : index = sample id, columns = [gene_loci_name], values = counts

        metaDF : columns = [sample_id, diagnosis]

        elifeLungGenes : list of ensembl gene names

        missingElifeGenes

        mapDF : columns = [HUGO_v35	ENSG_v35	ENSG_v39]

    '''
    logger.info( "BEGIN" )


    transposedCountsDF = loadCounts()
    logger.info( f'transposedCountsDF.shape : {transposedCountsDF.shape}' )
    logger.info( f'{transposedCountsDF.iloc[0:5, 0:5]}' )

    metaDF = loadMetaData()
    logger.info( f'metaDF.shape : {metaDF.shape}' )
    logger.info( f'metaDF.head :\n{metaDF.head()}' )

    elifeLungGenes, missingElifeGenes, mapDF = selectFeatures( HUGO_featuresNames )
    logger.info( f'len(elifeLungGenes) {len(elifeLungGenes)}' )
    logger.info( f'first 3 {elifeLungGenes[:3]} last 3 {elifeLungGenes[-3:]}' )

    logger.info ( f'missingElifeGenes : {missingElifeGenes}' )

    logger.info( f'mapDF.shape : {mapDF.shape}' )
    logger.info( f'mapDF.head():\n{mapDF.head()}' )

    # where did this Gene come from? Why is it missing
    logger.warning(f'AEDWIP where did this come form and why is it missing? ENSG00000182378.15_PAR_Y')
    biomakerGenes = [s for s in elifeLungGenes if s != "ENSG00000182378.15_PAR_Y"]
    XDF = transposedCountsDF.loc[:, biomakerGenes]
    logger.info( f'XDF.shape : {XDF.shape}' )
    logger.info( f'XDF.iloc[0:5, 0:3] : \n{XDF.iloc[0:5, 0:3]}' )

    logger.info( "END" )
    return (XDF, metaDF, elifeLungGenes, missingElifeGenes, mapDF)