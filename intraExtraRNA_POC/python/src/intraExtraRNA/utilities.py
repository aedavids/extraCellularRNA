'''
Created on May 31, 2023

@author: andrewdavidson
aedavids@ucsc.edu

public functions

def load(source, localCacheDir, verbose=False):

def selectSamples(colDataDF: pd.DataFrame, 
                  countDF:   pd.DataFrame, 
                  listOfCategories) -> pd.DataFrame:
                  
'''

import pathlib as pl
import pandas as pd
import shutil

################################################################################
def load(source, localCacheDir, verbose=False):
    '''
    ref: file:///Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/terra/jupyterNotebooks/cibersort/createCibersortMixtureMatrix.html

    reading large files over a NFS mount is slow. loadCache() will
    copy the source file into the local cache if it does not not already exist
    
    arguments:
        source:
            file path
        
        localCacheDir:
            path.
            
        verbose:
            if True will print the full local cache path to the file
    '''
    # we can not join, combine source if it start from the root of the file system
    tmpSource = source
    if source[0] == "/":
        tmpSource = source[1:]
    
    localTargetPath = pl.Path(localCacheDir,  tmpSource)
    if verbose:
        print("localTargetPath:\n{}\n".format(localTargetPath))
            
    localTargetPath.parent.mkdir(parents=True, exist_ok=True)

    if not localTargetPath.exists():
        if verbose:
            print("localTargetPath:{} does not exits".format(localTargetPath))
            print("copy {} to local cache".format(source))            
        #! cp $source $localTargetPath
        shutil.copy(source, localTargetPath)
        
    return localTargetPath

################################################################################
def selectSamples(colDataDF: pd.DataFrame, 
                  countDF:   pd.DataFrame, 
                  listOfCategories) -> pd.DataFrame:
    '''
    arguments:
        colDataDF:
            type pandas dataframe
            index is 'geneId'
                            
        countDF: 
            type pandas dataframe
            
            example groupedByGeneDF
            
        listOfCategories
            type list of strings
            
            example:
                ['PAAD', 'Pancreas']
                
    '''
    selectRows = colDataDF.loc[:, "category"].isin(listOfCategories)
    sampleIdsSeries = colDataDF.loc[selectRows, "sample_id"]
    
    # print(f'sampleIdsSeries.shape : {sampleIdsSeries.shape}')
    # display(sampleIdsSeries[0:5])
    
    sampleIdsList = sampleIdsSeries.to_list()
    #cols = ["geneId"] + sampleIdsList
    # print(f'cols[0:5] : {cols[0:5]}')
    
    # retDF = countDF.loc[:, sampleIdsList]
    retDF = countDF.loc[sampleIdsList, :]
    
    return retDF

################################################################################
def testLoad():
    source = "/private/groups/kimlab/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25/ciberSort/testSignatureGenes.txt"
    load( source )

################################################################################
if __name__ == '__main__':
    testLoad()
