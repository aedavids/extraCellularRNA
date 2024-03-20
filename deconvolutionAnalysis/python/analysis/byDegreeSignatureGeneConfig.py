#
# byDegreeSignatureGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

# ByDegreeSignatureGeneConfig command line display the doc string 
'''
The byDegreeSignatureGeneConfig class should not be used from cli.
Use the cli for testing varg parsing. see main()
'''

import ast
import logging
import os
import pandas as pd

from analysis.byDegreeSignatureGeneConfigCLI import ByDegreeSignatureGeneConfigCLI
from analysis.bestSignatureGeneConfig import BestSignatureGeneConfig
from analysis.utilities import fileNameToDictKey, findIntersectionsWithDegree, loadDictionary
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-12-17'
__updated__ = '2023-12-17'

##############################################################################
class ByDegreeSignatureGeneConfig(SignatureGeneConfiguration):
    '''
    genes that are shared between many categories do not make good discrimnators

    selects genes from intersections with degree = x

    ref: analysis.test.testByDegreeSignatureGeneConfig.py TestByDegreeSignatureGeneConfig
    '''

    logger = logging.getLogger(__name__)

    ################################################################################
    def __init__(self, 
                    dataSetName : str, 
                    design : str, 
                    padjThreshold : float, 
                    lfcThreshold : float, 
                    n : int, 
                    localCacheRootPath : str, 
                    title : str,
                    intersectionDictionaryPath : str,
                    degree : str,
                    )  :
        super().__init__(
                    dataSetName=dataSetName,
                    design=design,
                    padjThreshold=padjThreshold,
                    lfcThreshold=lfcThreshold,
                    n=n,
                    localCacheRootPath=localCacheRootPath,
                    title=title
        )
    
        '''
        see pipeline.dataFactory.signatureGeneConfig.SignatureGeneConfig().__init__ 
        for documentation

        additional arguments:
            intersectionDictionaryPath : str
                output from upstream stage of pipeilne. A dictionary key is multi-index identifying sets 
                that share elements. Value is the list of shared elements. 

                format example: ./analysis/test/data/intersection.dict

                see plots.test.testUpsetPlots testIntersection()

            degree : int
                see findGenes()
        '''
        self.logger.info(f'BEGIN')

        self.intersectionDictionaryPath = intersectionDictionaryPath
        self.degree = degree

        self.intersectionDict = loadDictionary(self.intersectionDictionaryPath)
        self.degreeDict = findIntersectionsWithDegree( self.intersectionDict, self.degree)


        # AEDWIP we are really only instrested in degree1
        # do not return all. makes upset plot uses less only one intersections with degree 83
        #
        # # create a set of all genes in interesections with degree
        # geneSet = set()
        # for key, genesLst in self.degreeDict.items():
        #     self.logger.debug(f'key: {key}, geneLst: {genesLst} geneSet: {geneSet}')
        #     s = set(genesLst)
        #     geneSet = geneSet.union( s )
        
        # self.genes = sorted( list(geneSet) )
        # self.logger.info(f'self.genes: {self.genes}')

        self.logger.info(f'END')

    ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        '''
        select genes in intersections with degree

        check log files to see which genes have been removed. They are probably 
        biologically interseting.
    
        arguments:
            deseqDF:
                results of DESeq2 as a pandas dataframe  

            fileName:
                useful for logging, debugging, and addition down stream processing   
    
        return:
            pandas dataframe
            
        ref: 
            BestSignatureGeneConfig.findGenes()

            findBestSignatureGenes(deseqDF, signatureGeneConfig)
                extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
        '''
        self.logger.info(f'BEGIN')
        self.logger.info(f'fileName: {fileName}')

        # do not select all cross all types. make upset plot useless
        # # retDF will contain all the unique genes discovered across all types
        # # so far
        # selectRows = deseqDF.loc[:, "name"].isin(self.genes)

        # default: return  an empty data frame
        retDF = pd.DataFrame(columns=['baseMean', 'lfcSE', 'log2FoldChange', 'name', 'padj', 'pvalue', 'stat'])
        category = fileName.removesuffix( "_vs_all.results" )
        key = fileNameToDictKey(fileName)
        
        if key in self.degreeDict: 
            elements = self.degreeDict[key]
            selectRows = deseqDF.loc[:, "name"].isin(elements)
            retDF = deseqDF.loc[selectRows, :]   
        
        numRows = retDF.shape[0]
        if retDF.empty or numRows == 0:
            self.logger.warning(f'fileName: {fileName} does not have any degree {self.degree} genes')     
#         bestDF = super().findGenes(deseqDF, fileName)

# # example of bestDF
# #          name     baseMean  log2FoldChange     lfcSE       stat         pvalue           padj
# # 6       C1orf61  2958.098269       -7.635107  0.380337 -20.074565   1.231553e-89   3.179694e-86
# # 3       ARPP21  1709.822182       -8.022019  0.399380 -20.086177   9.748529e-90   2.710541e-86
# # 8       OLIG1  1323.934363       -7.620726  0.421037 -18.099881   3.193580e-73   2.355819e-70

#         # select rows that are not poorDiscriminators
#         selectRows = ~bestDF.loc[:, "name"].isin(self.poorDiscriminators)
#         retDF = bestDF.loc[selectRows, :]

        self.logger.info(f'END')
        return retDF
    


################################################################################
def main(inCommandLineArgsList=None):
    '''
    use to test parsing of vargs
    '''
    # we only configure logging in main module
    loglevel = "WARN"
    loglevel = "INFO"
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    logger = logging.getLogger(__file__)
    logger.warning("BEGIN")

    #
    # always log run time env to make debug easier
    #
    ORIG_PYTHONPATH = os.environ['PYTHONPATH']
    logger.warning(f'PYTHONPATH : {ORIG_PYTHONPATH}')
    logger.warning(f'FILE: {__file__}')
    logger.warning(f'PWD: {os.getcwd()}')

    cli = ByDegreeSignatureGeneConfigCLI(version=__version__ , 
                            author=__author__ ,
                            date=__date__, 
                            update=__updated__)
    
    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    design          = cli.args.design 
    padjThreshold   = cli.args.padjThreshold
    lfcThreshold    = cli.args.lfcThreshold
    dataSetName     = cli.args.dataSetName
    number          = cli.args.number
    title           = cli.args.title
    localCacheRoot  = cli.args.localCacheRoot
    intersectionDictPath = cli.args.intersectionDict

    logger.warning(f'END')

################################################################################
if __name__ == '__main__':
    '''
    hack used to test parsing of vargs
    '''
    main()
    # main( inCommandLineArgsList=["--design", "tilda_gender_category",
    #                              "--padjThreshold", "0.001",
    #                              "--lfcThreshold", "2.0",
    #                              "--dataSetName", "GTEx_TCGA",
    #                              "--number" , "10",
    #                              "--title", "a vs all",
    #                              "--localCacheRoot", "myCacheRoot",
    #                              "--intersectionDict", "./analysis/test/data/intersection.dict",
    #                             ])


