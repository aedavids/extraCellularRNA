#
# plasmaGAN.py
# Andrew E. Davidson
# aedavids@ucsc.edu
#
'''
public functions 
   
'''

import logging
import pandas as pd

from intraExtraRNA.elifeUtilities import loadCounts
from intraExtraRNA.elifeUtilities import loadMetaData
from intraExtraRNA.elifeUtilities import selectFeatures

logger = logging.getLogger(__name__)

def loadCountData(
    HUGO_featuresNames : list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str], pd.DataFrame]:
    '''
    TODO

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
    biomakerGenes = [s for s in elifeLungGenes if s != 'ENSG00000182378.15_PAR_Y']
    XDF = transposedCountsDF.loc[:, biomakerGenes]
    logger.info( f'XDF.shape : {XDF.shape}' )
    logger.info( f'XDF.iloc[0:5, 0:3] : \n{XDF.iloc[0:5, 0:3]}' )

    return (XDF, metaDF, elifeLungGenes, missingElifeGenes, mapDF)