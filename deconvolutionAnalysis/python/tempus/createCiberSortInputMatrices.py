#
# createCiberSortInputMatrices.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref:
# extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
# extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/cibersortMixtureMatrixFactory.py
# # createCiberSortInputMatrices display the doc string 
'''
Functions to create input matrices for CiberSort from count data created by 
"Crate" , formally known as 'Complete Seq'
    
ref:
    extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
    extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/cibersortMixtureMatrixFactory.py
'''

import logging
logger = logging.getLogger(__file__)
import pandas as pd

def createSignatureMatrix(
        genesOfInterest : list[str],
        normalizedCountDF : pd.DataFrame,
        metaDF : pd.DataFrame,
        useMedian : bool = False
    ) :
    '''
    create a signature matrix for CiberSort
    
    arguments:
        genesOfInterest : list[str]

        normalizedCountDF : pd.DataFrame
            The DESeq2 normalized count data. Each row is a gene.
            the index column should be named 'gene_id'

            sample data:
                gene_id, gene_name, gene_biotype, SLDK3_T1_100_S1_L007
                (A)n, (A)n, Microsatellite, 0.5450089567650692
                (AAA)n, (AAA)n, Microsatellite, 0
                (AAAAAAC)n, (AAAAAAC)n, Microsatellite, 0
                (AAAAAAG)n, (AAAAAAG)n, Microsatellite, 0

        metaDF : pd.DataFrame,
            the DESeq2 colData. This data describes each sample. 
            E.G. sample_id category, ...
            The index column should be note be named 'sample_id' or 'category'

            ./illumina/20241107/raw/metadata.csv

            sample data:
                sample_id              category
                SLDK3_T1_100_S1_L007   A
                SLDK3_T1_100K_S2_L007  B
                SLDK3_T1_10K_S3_L007   A
                SLDK3_T1_1K_S4_L007    B

        useMedian : bool
            if True use the median to calculate the signature gene mean value for each category
            if False use the mean to calculate the signature gene mean value for each category
        
    returns:

    ref:
        extraCellularRNA/deconvolutionAnalysis/python/tempus/test/testCreateCiberSortInputMatrices.py
    '''
    logger.info("BEGIN")

    # sort the gene list and remove any duplicates
    sortedGeneList = sorted( set(genesOfInterest) ) 


    transposedDF = normalizedCountDF.transpose(copy=True)
    logger.info(f'transposeGroupByDF\n{transposedDF}')   

    # join the metaDF, we need the 'category' col so we can
    # calculate the  signature gene mean value for each category 
    joinDF =  pd.merge(left=transposedDF, 
                        right=metaDF.loc[:,["sample_id", "category"]], 
                        how='inner', 
                        left_index=True, 
                        right_on="sample_id")      
    logger.info(f'joinDF:\n{joinDF}')

    # calculate the expected values for each category  
    genesDF = joinDF.loc[ :, sortedGeneList + ["category"] ] 
    if useMedian :
        # weird duplicated log so I can set debugger break points
        logger.info(f'useMedian : {useMedian} calling median()')
        signatureDF = genesDF.groupby("category").median()
    else:
        logger.info(f'useMedian : {useMedian} calling mean()')
        signatureDF = genesDF.groupby("category").mean()
    
    # convert to cibersort expected upload format
    ciberSortSignatueDF = signatureDF.transpose(copy=True)
    ciberSortSignatueDF.index.name = "gene_id"
    
    # weird. cciberSortSignatueDF.columns.name = 'category'. This name is not
    # saved by pd.to_csv(). Set to none to make it easier to write unit test
    ciberSortSignatueDF.columns.name = None

    logger.info(f'END')
    return ciberSortSignatueDF