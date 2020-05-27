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
    classdocs
    '''
    
    logger = logging.getLogger(__name__)


    ################################################################################
    def __init__(self, inputPath, outputDirPath, adjPValue, logFoldChange):
        '''
        aedwip
        '''
        self.inputPath = Path(inputPath)
        self.outputDir = Path(outputDirPath)
        self.adjPValue = adjPValue
        self.logFoldChange = logFoldChange
        
        self.outputDir.mkdir(parents=True, exist_ok=True)
        
        inputFileName= PurePath(self.inputPath).name
        self.outputFile = self.outputDir.joinpath( inputFileName )
        
        # sample data
#         "","baseMean","log2FoldChange","lfcSE","pvalue","padj"
#         "7SK",37.988818348561,-0.203567415448379,0.412184624016707,0.386797918003106,0.736868517658304
#         "A1BG",57.8308345209208,-0.351497983353878,0.549275107282857,0.132784346152779,0.493600161162662
#         "A2M",8.95172936552604,-0.105027332164006,0.466985799757005,0.39752525626016,0.742161089145889

        LOG_IDX = 2
        P_ADJ_IDX = 5
        with open(self.inputPath, mode="r") as inputFH:
            with open(self.outputFile, mode="w+") as outputFH:
                headers = inputFH.readline()
                outputFH.write(headers) 
                self.logger.info("headers:{}".format(headers))
                
                for line in inputFH:
                    tokens = line.strip().split(',')
                    fancy = np.array(tokens)[ [LOG_IDX, P_ADJ_IDX] ]
                    if "NA" in fancy:
                        # array length is 2, 'in' will be fast 
                        continue
                    
                    logFold, padj =  fancy.astype(float)
                    if logFold >= self.logFoldChange and padj >= self.adjPValue:
                        outputFH.write(line)
                        self.logger.info("line:{}".format(line))
            outputFH.close()
        inputFH.close()
        
    ################################################################################
    def getOutputFileName(self):
        return self.outputFile