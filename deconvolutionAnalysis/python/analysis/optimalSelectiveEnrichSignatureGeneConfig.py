#
# optimalSelectiveEnrichSignatureGeneConfig.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

# OptimalSelectiveEnrichSignatureGeneConfigCLI CommandLine display the doc string 
'''
The OptimalSelectiveEnrichSignatureGeneConfig class should not be used from cli.
Use the cli for testing varg parsing. see main()
'''

import glob
import logging
import os
import pandas as pd
from plots.upsetPlots import UpsetPlot
import pprint as pp
import traceback

from analysis.bestSignatureGeneConfig import BestSignatureGeneConfig
from analysis.optimalSelectiveEnrichSignatureGeneConfigCLI import OptimalSelectiveEnrichSignatureGeneConfigCLI
from analysis.utilities import findAllGenes
from analysis.utilities import fileNameToDictKey
from analysis.utilities import findIntersectionsWithDegree
from analysis.utilities import loadDictionary
from analysis.utilities import loadSet
from pipeline.dataFactory.driver import _countExtraHeaderLines
from pipeline.dataFactory.driver import runSelectGenesOfInterest
from pipeline.dataFactory.upsetPlotDataFactory import UpsetPlotDataFactory

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2024-01-04'
__updated__ = '2024-01-04'

##############################################################################
class OptimalSelectiveEnrichSignatureGeneConfig(BestSignatureGeneConfig):
    '''
    genes that are shared between many categories do not make good discrimnators

    find best "n" genes. see BestSignatureGeneConfig.findGenes() for details
    then remove genes that are elements of intersections with high degree.
    The degree of a set is the number of set that have elements in the intersection

    ref: analysis/test/testOptimalSelectiveEnrichSignatureGeneConfig.py
    '''

    logger = logging.getLogger(__name__)

    ################################################################################
    def __init__(self, 
                    dataSetName : str, 
                    design : str, 
                    padjThreshold : float, 
                    lfcThreshold : float, 
                    windowLength : int, 
                    localCacheRootPath : str, 
                    title : str,
                    deseqResultsDir : str,
                    #upstreamDeseqResultsDir : str, do we need this? _step1 will pass this to findGenes()
                    categories : list[str],
                    startIdx : int,
                    maxNumberOfGenes : int,
                    upStreamIntersectionDictionaryPath : str,
                    historicGeneSetPath : str
                    )  :
        '''
        ref: analysis/test/testOptimalSelectiveEnrichSignatureGeneConfig.py

        Add degree1 genes to resuls from an upstream pipeline stage
        see pipeline.dataFactory.signatureGeneConfig.SignatureGeneConfig().__init__ 
        for documentation

        additional arguments:
            windowLength : int
                defines the range of rows to search for candidate genes
                The range will be startIdx ... startIdx + windowLength
                
            deseqResultsDir : str
                path to directory of deseq result files we will search for additional genes.
                This should probably be the original results file. i.e. not the output of previous 
                pipeline stages

            categories : list[str]
                list of classes, types, categories to enrich. If class does not exist in ??? historic ???
                It will enriched an added to output

            startIdx : int
                count start at zero
                startIdx can be used to reduce the search space

                once all candidate deseq results have been sorted. we start the search from this value
                The search will end after consdierting at most startIdx + n  rows.

                example if upstream stage was best100. startIdx should be 100

            maxNumberOfGenes : int
                for class in the list of categories
                    will try to add at most maxNumberOfGenes degree 1 genes
            
            upStreamIntersectionDictionaryPath : str
                output from upstream stage of pipeilne. Dictionary key is multi-index identifying sets that share elements.
                value is the list of shared elements

                used to find degree1 genes from previous upstream stages to add to

                example: ./analysis/test/data/intersection.dict

                see plots.test.testUpsetPlots testIntersection()

            historicGeneSetPath : str
                file containing a python set. use to ensure we do not add genes that had been considered or removed by 
                upstream stages

        '''        
        super().__init__(
                    dataSetName=dataSetName,
                    design=design,
                    padjThreshold=padjThreshold,
                    lfcThreshold=lfcThreshold,
                    n=windowLength,
                    localCacheRootPath=localCacheRootPath,
                    title=title
        )
        self.logger.info(f'BEGIN')

        self.deseqResultsDir            = deseqResultsDir
        # self.upstreamDeseqResultsDir    = upstreamDeseqResultsDir
        self.categories                 = categories

        self.startIdx                   = startIdx
        self.windowLength = windowLength
        self.endIdx = self.startIdx + self.windowLength #* len(self.candidateResultsFiles))
             
        self.maxNumberOfGenes           = maxNumberOfGenes
        self.upStreamIntersectionDictionaryPath = upStreamIntersectionDictionaryPath
        self.historicGeneSetPath = historicGeneSetPath

        self.upstreamIntersectionDict = loadDictionary(self.upStreamIntersectionDictionaryPath)
        # use upstream1Dict to find previouly discover degree1 genes
        self.upstreamDegree1Dict = findIntersectionsWithDegree(self.upstreamIntersectionDict, degree=1)

        # use historicGeneSet to ensure we do not add back genes
        # that had been previously removed or considered
        self.historicGeneSet = loadSet(self.historicGeneSetPath)

        # make sure all upstream streams are included in history
        # they should already be part of history. 
        # this just make it alittle easier to use this stage
        allUpstreamGenesSet = findAllGenes(self.upstreamIntersectionDict)
        self.historicGeneSet = self.historicGeneSet.union(allUpstreamGenesSet)

        # keep track of genes that have been added 
        self.locallyAddedGenes = []
     
        #
        # find degree1 gene candidate 
        #
           
        self.candidateResultsFiles  = self._findCandidateResultsFiles()
        self.logger.info(f'startIdx : {self.startIdx} n : {self.n} len(self.candidateResultsFiles) : {self.candidateResultsFiles}, endIdx : {self.endIdx}')
        self.candidateRankedDF = self._loadLocallyRankedResults()

        candidateIntersectionDict, candidateDegree1Dict = self._getCandidateDegree1Dict()
        self.candidateIntersectionDict = candidateIntersectionDict
        self.candidateDegree1Dict = candidateDegree1Dict

        self.logger.info(f'END')

    ################################################################################
    def findGenes(self, deseqDF : pd.DataFrame, fileName : str ) -> pd.DataFrame :
        '''
        TODO
        Find genes that that are statistically signifigant with  lfc <= threshold or >= threshold
        then removes genes from sets with high degree. These genes shared between many 
        categories and are not good descriminators.

        check log files to see which genes have been removed. They are probably 
        biologically interseting.
    
        arguments:
            deseqDF:
                results of DESeq2 as a pandas dataframe  
                This should be the output of an upstream pipeline stage

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
        
        # deseqDF is passed to use from upstreamPipeline _step1()
        # it is the output from an upstream stage
        # passing the original "pre pipeline " results file is a bug.
        # this class can not be a head stage

        category = fileName.removesuffix( "_vs_all.results" )
        key = fileNameToDictKey(fileName)
        

        # default:  pass through output from upstream stage
        retDF = deseqDF

        if category in self.categories:
            try:
                if not key in self.candidateDegree1Dict:
                    # return pass through
                    self.logger.warning(f'unable to enrich: no candidate degree 1 genes for {category}')
                    return retDF
                
                candidateGenes = self.candidateDegree1Dict[key] 
                newCandidateGenes = []

                # find genes that have not been considered before
                for newGene in candidateGenes:
                    if not newGene in self.locallyAddedGenes :
                        if not newGene in self.historicGeneSet:
                            newCandidateGenes.append( newGene )
                            self.locallyAddedGenes.append( newGene )
                self.logger.info(f'category: {category} newCandidateGenes : {newCandidateGenes}')

                # get the deseq result rows for the new candidates
                selectRows = self.candidateRankedDF.loc[:, "name"].isin(newCandidateGenes)
                candidateDF = self.candidateRankedDF.loc[selectRows, :]
                self.logger.warning(f'AEDWIP candidateDF:\n{candidateDF}')

                selectRows = candidateDF.loc[:, 'category'] == key
                candidateDF = candidateDF.loc[selectRows, :].sort_values(by="baseMean", ascending=False)
                self.logger.warning(f'AEDWIP sorted by baseMean candidateDF:\n{candidateDF}')

                if key in self.upstreamDegree1Dict :
                    # selectRows = deseqDF.loc[:, "name"].isin(newCandidateGenes)
                    # upstreamDegree1DF = deseqDF.loc[selectRows, :].sort_values(by="baseMean", ascending=False)
                    # nRows = upstreamDegree1DF.shape[0]
                    upstreamDegree1Genes = self.upstreamDegree1Dict[key]
                    nRows = len(upstreamDegree1Genes)
                    if  nRows < self.maxNumberOfGenes:
                        numAdd = self.maxNumberOfGenes - nRows
                        self.logger.warning(f'AEDWIP category : {category} self.maxNumberOfGenes : {self.maxNumberOfGenes} add {numAdd} genes')
                        addDeseqDF = candidateDF.head(n=numAdd)
                        self.logger.warning(f'AEDWIP addDeseqDF:]\n{addDeseqDF}' )
                        retDF = pd.concat( [retDF, addDeseqDF])
                else :
                    # upstream output did not have any degree1 genes for this category
                    addDeseqDF = candidateDF.head(n=self.maxNumberOfGenes)
                    retDF = pd.concat( [retDF, addDeseqDF] )
                    self.logger.warning(f'AEDWIP no upstream degree1 genes for {category} addDeseqDF:]\n{addDeseqDF}' )


            except Exception as e:
                self.logger.error(f'fileName : {fileName}')
                self.logger.error(f'exception : {e}')
                self.logger.error(traceback.print_exc() )

        # genes from self.candidateRankedDF have two columns we want to remove
        retDF = retDF.loc[:, ~retDF.columns.isin( ['categoryIdx', 'category'] ) ]
        self._checkForDuplicates(retDF)
        self.logger.warning(f'AEDWIP retDF:\n{retDF}')

        self.logger.info(f'END')
        return retDF

    ################################################################################
    def _checkForDuplicates(self, retDF : pd.DataFrame) :
        '''
        quality control. just print warning 
        '''
        self.logger.info('BEGIN')
        dupsSeries  = retDF.loc[:, ['name']].duplicated()
        if dupsSeries.sum() > 0:
            self.logger.error(f"duplicates found in retDF :\n{retDF.loc[dupsSeries, ['name']]}")

        self.logger.info('END')

    ################################################################################
    def _findCandidateResultsFiles(self) -> list[str]:
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
        for fPath in self.candidateResultsFiles:
            numRowsToSkip = _countExtraHeaderLines(fPath)
            self.logger.info(f"fPath: {fPath}")
            self.logger.info(f"numRowsToSkip: {numRowsToSkip}")
            deseqDF = pd.read_csv(fPath, skiprows=numRowsToSkip)

            fileName = os.path.basename(fPath)
            bioSignifigantDF = self._select(deseqDF, fileName)

            # make sure the bio sigfigant genes are in history
            # they should already be in history
            # self.logger.error("AEDWIP DO WE NEED TO DO THIS? fix unites test input")
            # nameIdx = bioSignifigantDF.columns.get_loc('name')
            # genes = bioSignifigantDF.iloc[0 : self.endIdx , nameIdx].to_list()
            # self.historicGeneSet = self.historicGeneSet.union( set(genes) )

            tmpDF = bioSignifigantDF.iloc[self.startIdx : self.endIdx , :]
            
            tmpDF = tmpDF.loc[:, ['name', 'baseMean', 'log2FoldChange', 'padj' ]]

            # add a category string constant to all rows
            # Whole_Blood_vs_all.results 
            #category = fileName.removesuffix( "_vs_all.results" )
            category = fileNameToDictKey(fileName)
            tmpDF['category'] = [category for i in range(tmpDF.shape[0])]

            retDF = pd.concat([retDF, tmpDF])

        retDF = retDF.reset_index(names="categoryIdx")
        self.logger.info(f'retDF.shape : {retDF.shape}')

        self.logger.info("END")
        return retDF

    ################################################################################   
    def _getBestSignatureGeneConfig(self) -> BestSignatureGeneConfig:
        '''
        broke out to make code more readable

        # we are derived from BestSignatureGeneConfig how ever we need to insure
        # BestSignatureGeneConfig.findGenes() is called        
        '''
        self.logger.info(f'BEGIN')

        sgc = BestSignatureGeneConfig(dataSetName=self.dataSetName,
                            design=self.design,
                            padjThreshold=self.padjThreshold,
                            lfcThreshold=self.lfcThreshold,
                            n=self.endIdx,
                            localCacheRootPath=self.localCacheRootPath,
                            title=self.title)
        
        self.logger.info(f'END')
        return sgc
    
    ################################################################################   
    def _getCandidateDegree1Dict(self) -> tuple[dict, dict]:
        '''
        broke out to improve read ablity

        ref: pipeline.upstreamPipeline
        '''
        self.logger.info(f'BEGIN')

        #
        # runSelectGenesOfInterest
        # filters/selects the genes of interest from each of the 1vsAll DESeq  results file
        # the selectedGeneSetsDict contains a key/value for each file
        # the outFileList is a file with a list of path to filtered DESeq results file.
        # ie only DESeq results we are interested
        #
     
        sgc = self._getBestSignatureGeneConfig()
        selectedGeneSetsDict, outFileList = runSelectGenesOfInterest(
                                            signatureGeneConfig=sgc, 
                                            candidateSignatureFileList=self.candidateResultsFiles,
                                            save=False,
                                            ) 
        
        upsetFactory = UpsetPlotDataFactory()
        upsetPlotDataDF, geneSetsDict = upsetFactory.createUpsetPlotDataFromSetDict(selectedGeneSetsDict)
        self.logger.error(f'AEDWIP upsetPlotDataDF:\n{pp.pformat(upsetPlotDataDF)}')

        upsetPlot = UpsetPlot(signatureGeneConfig=sgc, 
                        geneSetsDict=geneSetsDict, 
                        geneSetsUpsetPlotData=upsetPlotDataDF)
        
        candidateIntersectionDict = upsetPlot.findIntersectionElements()  
        candidateDegree1Dict = findIntersectionsWithDegree(candidateIntersectionDict, degree=1)

        self.logger.info(f'END')
        return (candidateIntersectionDict, candidateDegree1Dict)
    
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

    cli = OptimalSelectiveEnrichSignatureGeneConfigCLI(version="testVersion" , 
                            author="testAuthor",
                            date="testDate", 
                            update="testUpdate")
    
    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    design          = cli.args.design 
    logger.info(f'design : {design}')

    padjThreshold   = cli.args.padjThreshold
    logger.info(f'padjThreshold : {padjThreshold}')

    lfcThreshold    = cli.args.lfcThreshold
    logger.info(f'lfcThreshold : {lfcThreshold}')

    dataSetName     = cli.args.dataSetName
    logger.info(f'dataSetName : {dataSetName}')

    windowLength          = cli.args.windowLength
    logger.info(f'number : {windowLength}')

    title           = cli.args.title
    logger.info(f'title : {title}')

    localCacheRoot  = cli.args.localCacheRoot
    logger.info(f'localCacheRoot : {localCacheRoot}')

    # do not use variable name that are close to python key works
    categories       = cli.args.classes
    logger.info(f'categories : {categories}')

    historicGeneSetPath = cli.args.historicGeneSetPath
    logger.info(f'historicGeneSetPath : {historicGeneSetPath}')


    maxNumberOfGenes = cli.args.maxNumberOfGenes
    logger.info(f'maxNumberOfGenes : {maxNumberOfGenes}')

    resultsDir = cli.args.resultsDir
    logger.info(f'resultsDir : {resultsDir}')

    startIdx = cli.args.startIdx
    logger.info(f'startIdx : {startIdx}')

    upStreamIntersectionDictionaryPath = cli.args.upStreamIntersectionDictionaryPath
    logger.info(f'upStreamIntersectionDictionaryPath : {upStreamIntersectionDictionaryPath}')


    logger.warning(f'END')

################################################################################
if __name__ == '__main__':
    '''
    hack used to test parsing of vargs
    '''
    main( 
        inCommandLineArgsList=["--design", "tilda_gender_category",
                                 "--padjThreshold", "0.001",
                                 "--lfcThreshold", "2.0",
                                 "--dataSetName", "GTEx_TCGA",
                                 "--windowLength" , "10",
                                 "--title", "a vs all",
                                 "--localCacheRoot", "myCacheRoot",
                                "--classes", "classAAA classBBB",
                                "--historicGeneSetPath", "testHistoricGeneSetPath",
                                "--maxNumberOfGenes", "999",
                                "--resultsDir", "testresultsDir",
                                "--startIdx", "0",
                                "--upStreamIntersectionDictionaryPath", "testupStreamIntersectionDictionaryPath",
                                ]

            # inCommandLineArgsList=["--help"]                
        )

