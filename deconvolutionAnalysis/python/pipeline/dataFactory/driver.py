#
# driver.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# create to break circular import problem

import logging
import pandas as pd
import pathlib as pl
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

###############################################################################
def runSelectGenesOfInterest(
        signatureGeneConfig : SignatureGeneConfiguration, 
        candidateSignatureFileList : list[str],  
        ) -> tuple[dict, list[pl.PosixPath]]:
    '''
    finds candidate biomarker genes

    for each file in candidateSignatureFileList
        call SignatureGeneConfiguration.findGenes()
        save results to csv file in signatureGeneConfig.getLocalCachedDir()
    
    arguments:
        signatureGeneConfig : SignatureGeneConfiguration
            contains run parmeters
            implements def findGenes(self, deseqDF : pd.DataFrame ) -> pd.DataFrame :
            
        candidateSignatureFileList: 
            a list of file paths to candidate signature gene files filter/select           
            
    returns: (selectedDict, outFileList)
        selectedDict : dictionary
            key: candidate signature filename
            value: pandas dataframe with DESeq results for selected rows
                ie header columns name,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
        
        outFileList: list [pl.PosixPath]
            list of paths for saved file
    '''
    logger = logging.getLogger(__name__)
    retDict = {}
    retOutFileList = []
    outDir = pl.Path( signatureGeneConfig.getLocalCachedDir() )
    if not outDir.exists() :
        outDir.mkdir(parents=True, exist_ok=True)    

    for csgpFile in candidateSignatureFileList:
        numRowsToSkip = _countExtraHeaderLines(csgpFile)
        logger.info(f"runSelectGenesOfInterest() csgpFile: {csgpFile}")
        logger.info(f"runSelectGenesOfInterest() numRowsToSkip: {numRowsToSkip}")
        deseqDF = pd.read_csv(csgpFile, skiprows=numRowsToSkip)
        signatureGenesDF = signatureGeneConfig.findGenes(deseqDF)   

        fileName = csgpFile.split("/")[-1]

        outFilePath = outDir.joinpath( fileName )
        signatureGenesDF.to_csv(outFilePath, index=False)
        logger.warning("saved to file: {}".format(outFilePath))

        retDict[fileName] = signatureGenesDF
        retOutFileList.append(outFilePath)

    ret = (retDict, retOutFileList)
    return ret

###############################################################################
def _countExtraHeaderLines(file: str) -> int:
    '''
    Depending on how the DESeq results file was saved there can be a 
    variable number of header lines

    calculate the number of lines to skip
    '''
    ret = 0
    firstHeaderLine = "name,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj"
    with open(file) as fd:
        found = False
        while not found:
            line = fd.readline().strip()
            # some times the head column names are in quotes
            line = line.replace('"', "")

            if not line:
                found = True
            
            if line == firstHeaderLine:
                found = True
            else:
                ret += 1

    return ret


