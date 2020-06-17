'''
Created on May 26, 2020

@author: andrewdavidson
'''
import logging
import numpy as np
from pathlib import Path, PurePath

################################################################################
class DESeqSelect(object):
    '''
    public functions
        __init__(self, inputFilePath)
        saveRows(self, outputDirPath, adjPValue, logFoldChange)
    '''
    
    logger = logging.getLogger(__name__)
    GENE_NAME_IDX = 0
    BASE_MEAN_IDX = 1
    LOG_IDX = 2
    P_ADJ_IDX = 5

    ################################################################################
    def __init__(self, inputPath):
        '''
        aedwip
        '''
        self.inputPath = Path(inputPath)
        
    ################################################################################ 
    def readVolcanoPlotData(self):   
        '''
        returns: (geneNames, x, y)
            type: numpy array
                geneNames: list of strings
                baseMean
                x: log2(fold-change)
                y: -log10(p-value)
        '''
        geneNames = []
        
        # numpy append is much slower than python list append
        # convert to numpy before return
        baseMeanList = []
        xList = []
        yList = []

        with open(self.inputPath) as fd:
            hdr = fd.readline()
            for line in fd:
                tokens = line.strip().split(',')                
                fancy = np.array(tokens)[ [self.BASE_MEAN_IDX, self.LOG_IDX, self.P_ADJ_IDX] ]
                # array length is 2, 'in' will be fast                 
                if "NA" in fancy:
                    continue
                    
                baseMean, log2Fold, adjP =  fancy.astype(float)
                
                if adjP == 0.0:
                    # log(0) is undefined
                    adjP = 0.000001
                                    
                baseMeanList.append( baseMean )
                xList.append( log2Fold )
                yList.append( adjP )
                geneNames.append(tokens[self.GENE_NAME_IDX])

        baseMeanNp = np.array( baseMeanList )
        xNP = np.array( xList )
        yNP = np.log10( yList ) * -1.0 
        # convert to numpy so we can use fancy indexing
        geneNamesNP = np.array(geneNames)
        
        return (geneNamesNP,baseMeanNp, xNP, yNP)
        
    ################################################################################    
    def saveRows(self, outputDirPath, adjPValue, logFoldChange):
        '''
        save rows that have adjusted p-values >= adjPValue and logFoldChange >= logFoldChange
        
        returns:
            output file name. Type PathLib.Path
        sample data
        
        "","baseMean","log2FoldChange","lfcSE","pvalue","padj"
        "7SK",37.988818348561,-0.203567415448379,0.412184624016707,0.386797918003106,0.736868517658304
        "A1BG",57.8308345209208,-0.351497983353878,0.549275107282857,0.132784346152779,0.493600161162662
        "A2M",8.95172936552604,-0.105027332164006,0.466985799757005,0.39752525626016,0.742161089145889
        '''
        outputDir = Path(outputDirPath)
        
        inputFileName= PurePath(self.inputPath).name
        outputDir.mkdir(parents=True, exist_ok=True)
        outputFile = outputDir.joinpath( inputFileName )

        with open(self.inputPath, mode="r") as inputFH:
            with open(outputFile, mode="w+") as outputFH:
                headers = inputFH.readline()
                outputFH.write(headers) 
                self.logger.info("headers:{}".format(headers))
                
                for line in inputFH:
                    tokens = line.strip().split(',')
                    fancy = np.array(tokens)[ [self.LOG_IDX, self.P_ADJ_IDX] ]
                    if "NA" in fancy:
                        # array length is 2, 'in' will be fast 
                        continue
                    
                    logFold, padj =  fancy.astype(float)
                    if logFold >= logFoldChange and padj >= adjPValue:
                        outputFH.write(line)
                        self.logger.info("line:{}".format(line))
            outputFH.close()
        inputFH.close()
        
        return outputFile
        