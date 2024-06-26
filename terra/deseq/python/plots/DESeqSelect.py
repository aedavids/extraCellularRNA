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
    
    #"name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
    # lfcshrink does not have stat col. "name","baseMean","log2FoldChange","lfcSE",        "pvalue","padj"
    # ['"SFTPA1"',  name 0
    #  '3211.04578787063', baseMean 1
    #  '11.753362315757', log2 2
    #  '0.199827704361778', lfcse 3
    #  '58.8174815564017', stat 4
    #  '0', pvalue 5
    #  '0' padj] 6
    
    
   
    GENE_NAME_IDX = 0
    BASE_MEAN_IDX = 1
    LOG_IDX = 2
    STAT_IDX = 5
    P_ADJ_IDX = 6

    ################################################################################
    def __init__(self, inputPath):
        '''
        aedwip
        '''
        self.inputPath = Path(inputPath)
        
    ################################################################################ 
    def readVolcanoPlotData(self, numHeaderLines, hackPadjIndx=P_ADJ_IDX):   
        '''
        arguments:
            numHeaderLines: number of lines to skip before before data begins
            
            hackPadjIndx:
                default is 6
                DESeq lfcShrink padj column index should be 5. does not have stat column
            
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
            for i in range(numHeaderLines):
                hdr = fd.readline()
               
            # aedwip = 0 
            for line in fd:
                tokens = line.strip().split(',')
                
                # if aedwip < 10:
                #     print(tokens)
                #     aedwip += 1
                # else:
                #     return
                
                try :
                    fancy = np.array(tokens)[ [self.BASE_MEAN_IDX, self.LOG_IDX, hackPadjIndx] ]
                except BaseException as e:
                    self.logger.error("exc:{} tokens:{}".format(e, tokens))
                    continue
                
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
                # example gene name "GLCE"
                gn = tokens[self.GENE_NAME_IDX].strip()
                if gn.startswith('"'):
                    gn = gn[1:]
                if gn[-1] == '"':
                    gn = gn[:-1]
                geneNames.append( gn )

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
        
    ################################################################################    
    def loadDESeqResultsAsStrings(self, numHeaderLines):
        '''
        we use upsetPlot of find sets of signature genes that have over lapping genes
        Use deseqResults to load data that can be use to print of the DESeq stats for 
        genes in the intersection. We can use this information to hand tune our 
        signature Genes
        arguments:
            numHeaderLines: number of lines to skip before before data begins
            
        returns: set
                key = geneNames
                value = list of values: the corresponding data for the 
                        gene in a DESeq result file. Values are Strings 
                        representations of real numbers
        '''
        retSet = {}
        
        with open(self.inputPath) as fd:
            for i in range(numHeaderLines):
                hdr = fd.readline()
               
            for line in fd:
                tokens = line.strip().split(',')
                key = tokens[0]
                value = tokens[1:]
                retSet[key] = value
        
        return retSet       
