'''
Created on May 31, 2023

@author: andrewdavidson
aedavids@ucsc.edu
'''
import pandas as pd
import pathlib

################################################################################
# type hint for standard collection requires python 3.9 GeneList = list[str]
#LoadTx2GenesTuple = tuple(pd.DataFrame, GeneList)
# def loadTx2Genes(filePath : pathlib.Path) -> LoadTx2GenesTuple:
def loadTx2Genes(filePath : pathlib.Path) -> (pd.DataFrame, list):
    '''
    assume filePath file is a two column file, first column is a transcript id, the second is a gene id
    The top of the file contains gene code entries of the form
    
    ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|,DDX11L1
    
    The bottom of the file has entries created from the UCSC genome browser repeat mask annotations
    
    example
    hg38_rmsk_ALR/Alpha_range=chr22_KI270739v1_random:50071-73985_5'pad=0_3'pad=0_strand=-_repeatMasking=none,ALR/Alpha

    
    arguments
        filePath : pathlib.Path
            example 
            Path("/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv")
            
    returns:
        #LoadTx2GenesTuple :tuple(pd.DataFrame, GeneList)
        (pd.DataFrame, List, list)
        
        tx2GeneDF: type pandas Dataframe.
                    pandas read_csv(filePath, names=['uberId', 'geneId'], index_col='geneId')
             
        list : gene code gene names
            example:
                ['AC010319.1', 'AC009318.3', 'AL358334.1', 'SMIM34A', 'RPL3P1']
               
        list: TE 'gene names'
            these may include repeated sequences
            example:
                ['(TGTACATG)n', '(AGACAGA)n', 'AluYk11', '(CACAT)n', '(GGCTCGT)n']
    '''
    tx2GeneDF = pd.read_csv(filePath, names=['uberId', 'geneId'], index_col='geneId')
    selectGeneCode = tx2GeneDF.loc[:, 'uberId'].str.startswith('ENST')
    gencodeDF = tx2GeneDF[selectGeneCode]
    gencodeGeneIdsMasterList = gencodeDF.index.values
    #selectTEs = ~tx2GeneDF.loc[:, 'uberId'].str.startswith('ENST')
    #txTE2GeneDF = tx2GeneDF[selectTEs]
    txTE2GeneDF = tx2GeneDF[~selectGeneCode]
    TE_geneIdsMasterList = txTE2GeneDF.index.values
    
    uniqueGencodeList = list(set(gencodeGeneIdsMasterList))
    uniqueTeList = list(set(TE_geneIdsMasterList))
    
    return (tx2GeneDF, uniqueGencodeList, uniqueTeList)
    
#     # bestGenesList geneId HUGO format.
# # countDF gene ids are in in ENSG format or loci . ex. ENSG00000227232.5 and (AAAAAC)n
# file="/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv"
# tx2GeneDF = pd.read_csv(file, names=['uberId', 'geneId'], index_col='geneId')
#
# selectTEs = ~tx2GeneDF.loc[:, 'uberId'].str.startswith('ENST')
# txTE2GeneDF = tx2GeneDF[selectTEs]
# print(txTE2GeneDF.shape)
#
# display(txTE2GeneDF.head(n=3))
#
# print()
# display(txTE2GeneDF.tail(n=3))
#
# TE_geneIdsMasterList = txTE2GeneDF.index.values
# print(len(TE_geneIdsMasterList))
# print(TE_geneIdsMasterList[0:3])
# print(TE_geneIdsMasterList[-3:])

################################################################################
if __name__ == '__main__':
    pass