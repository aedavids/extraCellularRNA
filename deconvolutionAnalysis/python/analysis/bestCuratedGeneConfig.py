#
# bestCuratedGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

# BestCuratedGeneConfig.pyCLI CommandLine display the doc string 
'''
The BestCuratedGeneConfig.py class should not be used from cli.
Use the cli for testing varg parsing. see main()
'''

import logging
import os
import pandas as pd

from analysis.bestCuratedGeneConfigCLI import BestCuratedGeneConfigCommandLine 
#from analysis.bestSignatureGeneConfig import BestSignatureGeneConfig

from analysis.utilities import fileNameToDictKey, findIntersectionsWithDegree, loadDictionary
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2024-01-08'
__updated__ = '2024-01-08'

##############################################################################
# class BestCuratedGeneConfig(BestSignatureGeneConfig):
class BestCuratedGeneConfig(SignatureGeneConfiguration):
    '''
    select top n genes sorted by base mean from the degree1 intersections

    psudo code
    - load interesectionDictPath
    - degree1Dict = findIntersectionsWithDegree(interesectionDict, degree=1
    )
    - findGenes(deseqResults, fileName)
        if category in degree1Dict
            genes = degree1Dict[category]
            resultsDF = deseqResuts[genes]
            return resultsDF.sort_values(by'baseMean', accending=False).head(n)
        else 
            return empty dataframe

    Use Case 1: run our training pipeline on hand crafted gene signature matrix
        - interesection dictionary degree1 geneLists define the gene signature matrix
        we want to train with. It was edited by hand or create programmatically. 

        - set n to very large number like 5,000 to cause all the degree1 values to be returns

        - pipeline stage driver shell script config
            + deseqResultsDir contains the results files that will be based to findGenes()
        
            + example: extraCellularRNA/deconvolutionAnalysis/bin/1vsAll-~gender_category/best10CuratedDegree1.sh
                * rootDir="/private/groups/kimlab/GTEx_TCGA"
                * deseqResultsDir="${rootDir}/1vsAll"

    Use Case 2: Automatially add genes to a signature matrix defined by the interesectionDict
        - interesection dictionary degree1 geneLists define the gene signature matrix
        we want to train with. 

        - set n to the number of degree 1 values to return. Typically n is a small number

        example n = 3 will return at most 3 genes from the categories' degree 1 genes

        see pipeine stage dirver instructions above
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
                    interesectionDictPath : str,
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

        interesectionDictPath
        '''

        self.logger.info(f'BEGIN')

        self.interesectionDictPath = interesectionDictPath;
        self.logger.info(f'self.interesectionDictPath:\n{self.interesectionDictPath}')
        self.intersectionDict = loadDictionary( self.interesectionDictPath )
        self.degree1Dict = findIntersectionsWithDegree( self.intersectionDict, degree=1 )

        self.logger.info(f'END')

   ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        '''
        TODO
        Find genes that that are statistically signifigant with  lfc <= -2.0 or >= 2.0

        arguments:
            deseqDF:
                results of DESeq2 as a pandas dataframe            

                
            fileName:
                useful for logging, debugging, and addition down stream processing   

        return:
            pandas dataframe
            if fileName category is in degree1 intersection dict
                return top N biologically signifigant genes sorted by baseMean
            else
                return top N from deseq Argument
            
        ref: findBestSignatureGenes(deseqDF, signatureGeneConfig)
            extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
        '''
        self.logger.info("BEGIN")

        key = fileNameToDictKey(fileName)
        self.logger.debug(f'AEDWIP fileName : {fileName} key: {key}')
        self.logger.info(f'key: {self.degree1Dict.keys()}')

        # make sure values are sorted by baseMean
        deseqDF = deseqDF.sort_values(by="baseMean", ascending=False)

        retDF = None
        if key not in self.degree1Dict:
            #
            # method 1
            # return desqDF did not provide vest over all results
            #
            # depending on the source of the deseqDF, the rows may be in alphabetic order
            # self.logger.warn(f'no degree1 genes found for category {key} returning deseqDF.head(n={self.n})')
            # retDF = deseqDF.head(n=self.n)

            #
            # method 2: 
            # default: return  an empty data frame
            # do not artifically add genes to curated signature/intersection dictionary
            retDF = pd.DataFrame(columns=['baseMean', 'lfcSE', 'log2FoldChange', 'name', 'padj', 'pvalue', 'stat'])            
            retDF = pd.DataFrame(columns=deseqDF.columns)    
            self.logger.warning(f'key: {key} not found in degree1Dict. returning empty data frame') 

            #
            # method 3
            #
            # retDF = super().findGenes(deseqDF, fileName)        
            # sortedD1DF = d1DF.sort_values(by='baseMean')
            # retDF = sortedD1DF.head(n=self.n)

            #
            # best method commit id ce467ff
            # I would not think this works in general
            # topN=10
            # upstreamTopN=500
            # degree=1
            # windowLength=500
            # upstreamRun="best${upstreamTopN}FindAllDegree${degree}_wl${windowLength}"
            # upstream="/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/${upstreamRun}"
            # upsetOut="${upstream}/training/${upstreamRun}.sh.out/upsetPlot.out"
            # intersectionDict="${upsetOut}/best${upstreamTopN}_findAllDegree${degree}_wl${upstreamTopN}.intersection.dict"

            retDF = deseqDF.head(n=self.n)

        else :
            # the stage that created the degree1 intersection dictionary was probably 
            # best500FindAllDegree1_wl500Degree1. It filtered by lfc and padj and sorted by baseMean
            # It is possible that the degree1 genes are outside of the lfc and padj ranges
            # it possible baseMean is small
            degree1Genes = self.degree1Dict[key]
            selectedRows = deseqDF.loc[:, 'name'].isin(degree1Genes)
            d1DF = deseqDF.loc[selectedRows, :]

            # we want to select degree1 genes that have the strongest signal
            # if we think about how we fit models. large baseMean should make
            # linear regresion miss classifcation error residual greatly. This should improve model
            # results.
            sortedD1DF = d1DF.sort_values(by='baseMean', ascending=False)
            retDF = sortedD1DF.head(n=self.n)

                      

        self.logger.info("END")
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

    logger = logging.getLogger(os.path.basename(__file__))
    logger.warning("BEGIN")

    #
    # always log run time env to make debug easier
    #
    ORIG_PYTHONPATH = os.environ['PYTHONPATH']
    logger.warning(f'PYTHONPATH : {ORIG_PYTHONPATH}')
    logger.warning(f'FILE: {__file__}')
    logger.warning(f'PWD: {os.getcwd()}')

    cli = BestCuratedGeneConfigCommandLine(version=__version__ , 
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
    interesectionDictPath = cli.args.interesectionDictPath

################################################################################
if __name__ == '__main__':
    '''
    hack used to test parsing of vargs
    '''
    main( inCommandLineArgsList=["--design", "tilda_gender_category",
                                 "--padjThreshold", "0.001",
                                 "--lfcThreshold", "2.0",
                                 "--dataSetName", "GTEx_TCGA",
                                 "--number" , "10",
                                 "--title", "a vs all",
                                 "--localCacheRoot", "myCacheRoot",
                                 "--interesectionDictPath", "myInteresectionDictPath",
                                ])


