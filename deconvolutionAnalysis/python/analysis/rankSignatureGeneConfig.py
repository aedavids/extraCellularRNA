#
# RankSignatureGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

# RankSignatureGeneConfigCLI CommandLine display the doc string 
'''
AEDWIP
'''
import glob
import logging
import os
import pandas as pd
import pprint as pp

from analysis.rankSignatureGeneConfigCLI import RankSignatureGeneConfigCommandLine
from analysis.utilities import fileNameToDictKey
from pipeline.dataFactory.driver import _countExtraHeaderLines
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-12-28'
__updated__ = '2023-12-28'

##############################################################################
class RankSignatureGeneConfig(SignatureGeneConfiguration):
    '''
    overview:
    1. concat all deseq results files
    2. select rows that are biologically signigant
    3. sort by base mean
    4. rows with strongest biologically signigant signal will be on top
    5. itterate over the sorted list
        if we find a "new" gene assign it to category with ie highest base mean.
        i.e. signal strength

    Assigning genes to category based on signal strength may still produce poor
    discriminators. It is possible genes discovered this way are members of 
    intersections with high degree
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
                    deseqResultsDir : str
                    )  :
        '''
        overview:
        1. concat all deseq results files
        2. select rows that are biologically signigant
        3. sort by base mean
        4. rows with strongest biologically signigant signal will be on top
        5. itterate over the sorted list
            if we find a "new" gene assign it to category with ie highest base mean.
            i.e. signal strength

        Assigning genes to category based on signal strength may still produce poor
        discriminators. It is possible genes discovered this way are members of 
        intersections with high degree
        '''
        super().__init__(
                        dataSetName=dataSetName,
                        design=design,
                        padjThreshold=padjThreshold,
                        lfcThreshold=lfcThreshold,
                        n=n,
                        localCacheRootPath=localCacheRootPath,
                        title=title
            )
        self.logger.info(f'BEGIN')

        self.deseqResultsDir = deseqResultsDir

        # list of all selected genes
        self.rankedGeneSet = set()

        self.resultsFiles  = self._findResultsFiles()
        self.localRankedDF = self._loadLocallyRankedResults()

        # top rows have the strongest biologically signifigant signals
        self.rankedDF      = self.localRankedDF.sort_values(by="baseMean", ascending=False)
        self.rankedDict    = self._rank()

        self.logger.info(f'END')

   ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        '''
        TODO AEDWIP

        arguments:
            deseqDF:
                results of DESeq2 as a pandas dataframe            

                
            fileName:
                useful for logging, debugging, and addition down stream processing   

        return:
            pandas dataframe
            
        ref: findBestSignatureGenes(deseqDF, signatureGeneConfig)
            extraCellularRNA/terra/jupyterNotebooks/signatureGenesUpsetPlots.ipynb
        '''
        self.logger.info("BEGIN")

        category = fileName.removesuffix( "_vs_all.results" )
        key = fileNameToDictKey(fileName)
        names = self.rankedDict[key]
        selectRows = deseqDF.loc[:, "name"].isin(names)
        retDF = deseqDF.loc[selectRows, :]

        self.logger.info("END")
        return retDF

    #
    # ######### private functions ############
    #

    ################################################################################
    def _findResultsFiles(self) -> list[str]:
        '''
        # get a list of deseq results file we want to select gene signatures from
        '''
        self.logger.info("BEGIN")
        
        deseqResultFiles = os.path.join(self.deseqResultsDir, "*.results")
        deseqResultFileList = glob.glob(deseqResultFiles)
        self.logger.info(f'deseqResultFileList : \n{pp.pformat(deseqResultFileList)}')

        self.logger.info("END")
        return deseqResultFileList
    
    ################################################################################
    def _loadLocallyRankedResults(self) -> pd.DataFrame:
        '''
        load the deseq results files
        filter and sort
        append the top n to return
        '''
        self.logger.info("BEGIN")

        retDF = pd.DataFrame()
        for fPath in self.resultsFiles:
            numRowsToSkip = _countExtraHeaderLines(fPath)
            self.logger.info(f"fPath: {fPath}")
            self.logger.info(f"numRowsToSkip: {numRowsToSkip}")
            deseqDF = pd.read_csv(fPath, skiprows=numRowsToSkip)

            fileName = os.path.basename(fPath)
            bioSignifigantDF = self._select(deseqDF, fileName)

            # 20 is a fudge factor.
            # we want to make sure we have more than n genes to search
            # when n == 20 and fudge was == 2 only one class would up with 20 genes
            tmpDF = bioSignifigantDF.head(n=self.n * 20)
            
            tmpDF = tmpDF.loc[:, ['name', 'baseMean', 'log2FoldChange', 'padj' ]]

            # add a category string constant to all rows
            # Whole_Blood_vs_all.results 
            #category = fileName.removesuffix( "_vs_all.results" )
            category = fileNameToDictKey( fileName )
            tmpDF['category'] = [category for i in range(tmpDF.shape[0])]

            retDF = pd.concat([retDF, tmpDF])

        retDF = retDF.reset_index(names="categoryIdx")
        self.logger.info(f'retDF.shape : {retDF.shape}')

        self.logger.info("END")
        return retDF

    ################################################################################
    def _rank(self) -> dict[str, list[str]] :
        '''
        TODO AEDWIP
        '''
        self.logger.info(f'BEGIN')

        retDict = dict()

        nTotalGenes = len(self.resultsFiles) * self.n
        nRows = self.rankedDF.shape[0]
        self.logger.info(f'num rows to search : {nRows} will try to find nTotalGenes : {nTotalGenes}')
        count = 0
        for i in range(nRows):
            candidateRowSeries = self.rankedDF.iloc[i, :]
            candidateName = candidateRowSeries.loc["name"]

            if not candidateName in self.rankedGeneSet :
                category = candidateRowSeries['category']
                # self.rankedGeneSet.add(candidateName)
                if not category in retDict:
                    retDict[category] = [candidateName]
                    self.rankedGeneSet.add(candidateName)
                    count = count + 1

                elif len(retDict[category]) < self.n:
                    retDict[category].append(candidateName)
                    self.rankedGeneSet.add(candidateName)
                    count = count + 1

            if count > nTotalGenes:
                self.logger.info(f'i : {i} expected nTotalGenes : {nTotalGenes}')
                break

        if count < nTotalGenes:
            self.logger.warning(f'i : {i} found count : {count} < expected nTotalGenes : {nTotalGenes}')

        debugKeys = retDict.keys()
        self.logger.info(f'len(retDict.keys()): {len(debugKeys)}')
        for k in debugKeys:
            self.logger.info(f'key: {k} numGenes: { len(retDict[k]) }')

        self.logger.info(f'END')
        return retDict

    
    ################################################################################
    def _select(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        '''
        select genes by padj and logfold changes and order by baseMean
        '''
        self.logger.info(f'BEGIN')
        
        # colsToReturn = deseqDF.columns
    
        #
        # find statistically signifigant genes
        #
        selectSignificantRowsPS = deseqDF.loc[:,"padj"] < self.padjThreshold
        self.logger.info("number of genes with padj < {} : {}".format(self.padjThreshold,
                                                        selectSignificantRowsPS.sum()))

        #
        # use absolute value of log fold change to select best 
        # biologically signifigant genes
        # 
        significantDF = deseqDF.loc[selectSignificantRowsPS, :]
        absPS = significantDF['log2FoldChange'].abs()
        significantDF2 = significantDF.assign(absLog2FoldChange=absPS)    

        selectBestUpRegulatedRows = significantDF2.loc[:, 'absLog2FoldChange'] >= self.lfcThreshold
        significantDF3 = significantDF2.loc[selectBestUpRegulatedRows, :]
            
        significantDF3 = significantDF3.sort_values( by = ['baseMean'], ascending=False)

        self.logger.info(f'fileName : {fileName} retDF.shape : {significantDF3.shape} ')

        self.logger.info(f'END')
        return significantDF3
    
 
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

    cli = RankSignatureGeneConfigCommandLine(version=__version__ , 
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
    deseqResultsDir = cli.args.deseqResultsDir

    logger.info(f'END')
 


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
                                 "--deseqResultsDir", "myDeseqResultsDirPath"
                                ])


