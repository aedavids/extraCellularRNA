#
# cibersortSignatureFactory.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
#
# CibersortSignatureMatrixFactoryCLI display the doc string 
'''
    ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb

    Given a set of 1vsAll results, creates a scaled GeneSignatureProfile Matrix.

    all genes in the results file will be include in the GeneSignatureProfile

    See extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/test/data/ 
    for example input data file formats
'''

'''
aedwip fix doc string Converts a transcript aggregated by gene loci count
s to ciber sort mixture matrix format. See extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/test/data/ for example input data file formats
'''

import logging
import os
import pandas as pd
from pipeline.dataFactory.utilities import loadCache
from pipeline.dataFactory.cibersortSignatureMatrixFactoryCLI import CibersortSignatureMatrixFactoryCLI
import sys

__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-10-04'
__updated__ = '2023-10-04'

###############################################################################
class CibersortSignatureMatrixFactory( object ):
    '''
    ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb

    Given a set of 1vsAll results, creates a scaled GeneSignatureProfile Matrix.

    all genes in the results file will be include in the GeneSignatureProfile

    public functions
        __init__
        getCiberSortSignatueDF
        save
    '''
    logger = logging.getLogger(__name__)

    ################################################################################    
    def __init__(self, 
                geneSignatureProfilesDataRootDir : str, 
                #oneVsAllDataDirDeprecated : str, 
                groupByGeneCountFilePath : str,   
                colDataFilePath : str,
                estimatedScalingFactorsFilePath : str,
                localCacheDir : str, 
                outdir  : str ="ciberSort",
                testSize : bool=None,
                verbose: bool=False,
                useMedian : bool=False) :
        '''
        arguments
            geneSignatureProfilesDataRootDir: str
                the path the the geneSignatureProfiles files. There should be file for each category.
                The file contains the DESeq2 results for the genes of interest.

                ex. "/GTEx_TCGA/geneSignatureProfiles/best/GTEx_TCGA_1vsAll-design:~__gender_+_category-padj:0.001-lfc:2.0-n:25"
                ```
                 $ wc -l Whole_Blood_vs_all.results
                26 Whole_Blood_vs_all.results
                
                $ head -n 10 Whole_Blood_vs_all.results
                name,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
                HBA2,939295.036780794,11.3962888262948,0.106055647081676,107.455747429632,0.0,0.0
                HBB,985081.569770207,11.3779865954891,0.11089940519715,102.597363577037,0.0,0.0
                HBA1,384428.792719921,11.3333149379708,0.107194068493031,105.727071444328,0.0,0.0
                ALAS2,3853.24560048757,10.9663473238829,0.122522623050016,89.5046731035637,0.0,0.0
                AC104389.6,1032.06794351778,10.6506713620897,0.215678197674272,49.3822346298297,0.0,0.0
                HBM,373.003020915427,10.6290392458312,0.152953063971155,69.492162954223,0.0,0.0
                HBD,1044.46629161082,10.1793798559813,0.140411756797509,72.4966348128612,0.0,0.0
                CXCR1,3110.51569343056,10.0579110598547,0.110749735982355,90.8165691831285,0.0,0.0
                HBG2,1115.92293113395,10.0184837892341,0.127200952521455,78.7610752171392,0.0,0.0
                ```

            # oneVsAllDataDirDeprecated: str
            #     the path to the results from 1vsAll results
                
            groupByGeneCountFilePath: str
                example: /private/groups/kimlab/GTEx_TCGA/groupbyGeneTrainingSets/GTEx_TCGA_TrainGroupby.csv
                ```
                $head -n 5 GTEx_TCGA_TrainGroupby.csv | cut -d , -f 1,2,3
                geneId,GTEX-1117F-0226-SM-5GZZ7,GTEX-1117F-0526-SM-5EGHJ
                (A)n,9,1
                (AAA)n,0,0
                (AAAAAAC)n,0,0
                (AAAAAAG)n,0,0
                ```
                
            colDataFilePath: str
                sample meta data
                example : /private/groups/kimlab/GTEx_TCGA/groupbyGeneTrainingSets/GTEx_TCGA_TrainColData.csv
                ```
                $ head -n 3 GTEx_TCGA_TrainColData.csv
                sample_id,participant_id,category,gender,age,dataSet
                GTEX-1117F-0226-SM-5GZZ7,GTEX-1117F,Adipose_Subcutaneous,Female,66.0,GTEx
                GTEX-1117F-0526-SM-5EGHJ,GTEX-1117F,Artery_Tibial,Female,66.0,GTEx
                ``` 
                
            estimatedScalingFactorsFilePath:
                DESeq normalizing factors. Scaling factor for each sample. Adjust for different
                library sizes and composition

                !!!! this file does not contain sampleId information.
                !!!! it is assumed to be in the same order as the samples in the groupByGeneCountFile

                ex. GTEx_TCGA/1vsAll/estimatedSizeFactors.csv

                ```
                $ head estimatedSizeFactors.csv 
                "sizeFactors(dds)"
                0.826165555539692
                0.683190314596956
                0.532358863795858
                0.837840308886721
                0.663139416826015
                ```
            
            localCacheDir : str
                path to local file cache
                example: /scratch/aedavids/tmp

                when reading files will first check localCache. If not in cache will copy
                ensures file are localized.

            outdir:
                string
                default = "ciberSort"
                    save() path will be "geneSignatureProfilesDataRootDir + "/" + outDir"
            testSize:
                integer
                default: None (i.e. select all)
                use to select a sub set for testing purposes
                
            verbose:
                boolean, default = False
                argument passed to loadCache()

            median:
                boolean, default = false
                if true use mean, else use mean to calculate expected count values for each gene in a category
        '''
        self.localCacheDir = localCacheDir

        self.geneSignatureProfilesDataRootDir = geneSignatureProfilesDataRootDir
        self.logger.info("geneSignatureProfilesDataRootDir\n{}".format(geneSignatureProfilesDataRootDir))
        
        # self.oneVsAllDataDirDeprecated = oneVsAllDataDirDeprecated
        # self.logger.info("oneVsAllDataDirDeprecated\n{}".format(oneVsAllDataDirDeprecated))
        
        self.groupByGeneCountFilePath = groupByGeneCountFilePath
        self.logger.info("groupByGeneCountFilePath\n{}".format(groupByGeneCountFilePath))

        self.outdir = outdir
        self.testSize = testSize
        self.verbose = verbose
        self.useMedian = useMedian
        
        # get a list of results files we want to use
        # and their deconvolution types. 
        self.suffix = "_vs_all.results"
        self.listOfResultsFiles = None
        self.listOfTypes = None
        self._initListOfTypes()
        
        # read the results into a dictionary of data frames
        self.resultsDFDict = None
        self._LoadResultsDFDict()
        
        # create a sorted list of all the uniqe signatue genes
        self.signatureGeneSet = None
        self._createsignatureGeneList()
        
        # free up memory
        self.resultsDFDict = None

        self.estimatedScalingFactorsFilePath = estimatedScalingFactorsFilePath;
        self.scalingFactors = None
        self._loadScalingFactors()     

        self.colDataFilePath = colDataFilePath
        self.colDataDF = None
        self._readColData()           
        
        self.groupByCountDF = None
        self._readGroupByCountDF()
        #self._readAndFilterGroupByCountDF()
        
        # self.colDataFilePath = colDataFilePath
        # self.colDataDF = None
        # self._readColData()
        
        # self.estimatedScalingFactorsFilePath = estimatedScalingFactorsFilePath;
        # self.scalingFactors = None
        # self._loadScalingFactors()
        
        self.ciberSortSignatueDF = None
        self._createSignatureMatrix()
        
        
        
    ################################################################################
    def getCiberSortSignatueDF(self):
        self.logger.info("BEGIN")
        return self.ciberSortSignatueDF
    
    def save(self, prefixStr=None) -> str:
        '''
        saves in a cibersort's expected format

        returns : str
            path to signature DF
        '''
        tOutDir = self.geneSignatureProfilesDataRootDir + "/" + self.outdir
        os.makedirs(tOutDir, exist_ok = True)
        
        if prefixStr:
            path = tOutDir + "/" + prefixStr + "SignatureGenes.tsv"
        else:
            path = tOutDir + "/signatureGenes.tsv"
        
        # self.ciberSortSignatueDF.to_csv(path, index=False, sep="\t")
        self.ciberSortSignatueDF.to_csv(path, sep="\t")
        self.logger.info("\nsaved to: {}".format(path))
        return path
 
        
    ################################################################################
    def _initListOfTypes(self):
        '''
        create a list of all the signature results file
        output by SignatureGenesUpsetPlots.ipynb
        '''
        self.logger.info("BEGIN")

        #listOfResultsFiles = !ls $bestDataRootDir
        #listOfResultsFiles = !ls $self.geneSignatureProfilesDataRootDir
        listOfResultsFiles = os.listdir( self.geneSignatureProfilesDataRootDir)
        if self.outdir in listOfResultsFiles:
            listOfResultsFiles.remove( self.outdir )

        # clean up listOfResults. all DESeq results file end with suffix '.results'
        tmpList = []
        for fileStr in listOfResultsFiles:
            if fileStr.endswith('.results'):
                tmpList.append( fileStr )

        listOfResultsFiles = tmpList

        self.listOfResultsFiles = sorted(listOfResultsFiles[0:self.testSize])
        self.logger.info("listOfResultsFiles:\n\t{}".format(self.listOfResultsFiles))

        # create list of types the we want to deconvolve the mixture matrix into
        self.listOfTypes = [ f.split(self.suffix)[0] for f in self.listOfResultsFiles]
        self.logger.info("listOfTypes:\n\t{}".format(self.listOfTypes))

    ################################################################################
    def _LoadResultsDFDict(self):
        '''
        read the DESeq 1vsAll results files into a dictionary of pandas data frames. 
        with 'type' string as key
        '''
        self.logger.info("BEGIN")

        self.resultsDFDict = dict()
        for i in range(len(self.listOfResultsFiles)):
            resultFile = self.listOfResultsFiles[i]
            deconvolutionType = self.listOfTypes[i]
            path = self.geneSignatureProfilesDataRootDir + "/" + resultFile
            path = loadCache(path, self.localCacheDir, verbose=self.verbose)
            df = pd.read_csv(path, sep=",")
            self.resultsDFDict[deconvolutionType] = df
            
            
    ################################################################################    
    def _createsignatureGeneList(self):
        '''
        The count file may have genes are are not interest. Use this list to filter.
        Come types may share signature genes
        '''
        self.logger.info("BEGIN")
       
        signatureGeneSet = set()
        for deconvolutionType,df in self.resultsDFDict.items():
            name = df.loc[:,'name']
            signatureGeneSet.update(name)
        
        # keep in sort order. makes debug easier
        self.geneListsorted = sorted( list(signatureGeneSet) )[0:self.testSize]
        self.logger.info("\nnumber of signature genes:{}".format(len(self.geneListsorted)))
        #self.logger.info(self.geneListsorted)    

   ################################################################################        
    def _readGroupByCountDF(self):
        '''
        load the groupByGene count data
        '''
        self.logger.info("BEGIN")
        
        path = loadCache(self.groupByGeneCountFilePath, self.localCacheDir, verbose=self.verbose)
        self.logger.info(f"loadCache completed: path: {path}")

        df = pd.read_csv(path, sep=",")
        
        # set index to geneId. will make join easier' When we transpose
        # the data frame the index will become the column names
        df = df.set_index('geneId')

        self.logger.info("_readAndFilterGroupByCountDF() shape:{}".format(df.shape))
        self.logger.debug(df.iloc[0:3, 0:4])

        #
        # Do not sort!!!!!
        # The DESeq estimatedScalingFactor file does not contain sampleId information
        # It is assumed that the estimatedScalingFactor are in the same order as the countFile sampleIds 
        # 

        # # geneId is an indx
        # #df = df[ df.loc[:, "geneId"].isin(self.geneListsorted)]
        # df = df[ df.index.isin(self.geneListsorted)]
        
        #  # sort makes debug easier
        # #df = df.sort_values( by=["geneId"] )
        # df = df.sort_index(ascending=True)

        # self.logger.info("\nshape:{}\niloc[0:3, 0:4]")
        # self.logger.info(df.iloc[0:3, 0:4])

        self.groupByCountDF = df        
        
    ################################################################################        
    def _readColData(self):
        '''
        load the colData. We only need the sample_id and category columns 
        '''
        self.logger.info("BEGIN")

        path = loadCache(self.colDataFilePath, self.localCacheDir, verbose=self.verbose)
        self.colDataDF = pd.read_csv(path, sep=",").loc[:, ['sample_id', 'category']]
        self.logger.info("_readColData() self.colDataDF.shape:{}".format(self.colDataDF.shape))
        self.logger.info("debug[0:3, :]")
        #self.logger.debug( self.colDataDF[0:3, :] )
              

    ################################################################################                   
    def _loadScalingFactors(self):
        '''
        TODO
        '''
        self.logger.info("BEGIN")

        path = loadCache(self.estimatedScalingFactorsFilePath, self.localCacheDir, verbose=self.verbose)
        self.scalingFactors = pd.read_csv(path, sep=",")
        self.logger.info("_loadScalingFactors() self.scalingFactors.shape:{}".format(self.scalingFactors.shape))
        self.logger.debug("iloc[0:3, :]")
        self.logger.debug( self.scalingFactors.iloc[0:3, :] )
        
    ################################################################################                   
    # def _printInfo(self, str, df):
    #     self.logger.info("{}.shape:{}".format(str, df.shape))
    #     self.logger.info("{}.iloc[0:3, 0:4]".format(str))
    #     self.logger.info(df.iloc[0:3, 0:4])

################################################################################                   
    def _createSignatureMatrix(self):
        '''
        todo
        '''
        self.logger.info("BEGIN")

        # copy so we do not accidently change original groupDF
        self.logger.debug(f'self.groupByCountDF\n{self.groupByCountDF}')
        transposeGroupByDF = self.groupByCountDF.transpose(copy=True)
        #self.logger.info("_createSignatureMatrix()")
        # self._printInfo('transposeGroupByDF', transposeGroupByDF)
        self.logger.debug(f'transposeGroupByDF\n{transposeGroupByDF}')

        
        # normalize counts
        # element wise multiplication . use values to to multiply a vector
        #self._printInfo('scalingFactors', self.scalingFactors)
        self.logger.debug(f'self.scalingFactors\n{self.scalingFactors}')

        #self.logger.info("************** AEDWIP_transposeGroupByDF.csv")
        #transposeGroupByDF.to_csv("AEDWIP_transposeGroupByDF.csv")

        normalizedDF = transposeGroupByDF *  self.scalingFactors.values
        #self._printInfo('normalizedDF', normalizedDF)
        self.logger.debug(f'normalizedDF\n{normalizedDF}')


        # select genes of interest
        normalizedDF = normalizedDF.loc[ :,self.geneListsorted]
        #self._printInfo('normalizedDF', normalizedDF)
        self.logger.debug(f'normalizedDF\n{normalizedDF}')


        # join the colData, we need the 'category' col so we can
        # calculate the  signature gene mean value for each category 
        joinDF =  pd.merge(left=normalizedDF, 
                            right=self.colDataDF.loc[:,["sample_id", "category"]], 
                            how='inner', 
                            left_index=True, 
                            right_on="sample_id")      
        self.logger.debug(f'joinDF:\n{joinDF}')


        # calculate the expected values for each category      
        genesDF = joinDF.loc[ :,self.geneListsorted + ["category"] ] 
        if self.useMedian :
            # weird duplicated log so I can set debugger break points
            self.logger.info(f'useMedian : {self.useMedian} calling median()')
            signatureDF = genesDF.groupby("category").median()
        else:
            self.logger.info(f'useMedian : {self.useMedian} calling mean()')
            signatureDF = genesDF.groupby("category").mean()
        
        # convert to cibersort expected upload format
        ciberSortSignatueDF = signatureDF.transpose(copy=True)
        ciberSortSignatueDF.index.name = "name"
        
        self.ciberSortSignatueDF = ciberSortSignatueDF
        # weird. cciberSortSignatueDF.columns.name = 'category'. This name is not
        # saved by pd.to_csv(). Set to none to make it easier to write unit test
        self.ciberSortSignatueDF.columns.name = None

        self.logger.info("END")        

################################################################################
def testLogCLI(logger, cli):
    logger.info("BEGIN")
    geneSignatureProfilesDataRootDir = cli.args.geneSignatureProfilesDataRootDir
    logger.info(f'geneSignatureProfilesDataRootDir : {geneSignatureProfilesDataRootDir}')

    groupByGeneCountFilePath = cli.args.groupByGeneCountFilePath
    logger.info(f' groupByGeneCountFilePath : {groupByGeneCountFilePath}')

    colDataFilePath = cli.args.colDataFilePath
    logger.info(f' colDataFilePath : {colDataFilePath}')

    scalingFactorsPath = cli.args.scalingFactorsPath
    logger.info(f' scalingFactorsPath : {scalingFactorsPath}')

    localCacheDir = cli.args.localCacheDir
    logger.info(f' localCacheDir : {localCacheDir}')

    outDir = cli.args.outDir
    logger.info(f' outDir : {outDir}')

    if cli.args.useMedian:
        useMedian = True
    else :
        useMedian = False
    logger.info(f' useMedian : {useMedian}')

    logger.info("END")

################################################################################
def testCLI(logger):
    logger.info("BEGIN")

    cli = CibersortSignatureMatrixFactoryCLI( 
                                             version=__version__ , 
                                             author=__author__ ,
                                             date=__date__, 
                                             update=__updated__)

    vargs= ["-varg1", "v1ArgData" ,
            "-varg2", "v2ArgData" ]

    # test backwards compatiblity, useMedian is not passed
    inCommandLineArgsListMean = [
                    "--groupByGeneCountFilePath", "groupByGeneCountFilePathArg",
                    "--geneSignatureProfilesDataRootDir", "geneSignatureProfilesDataRootDirArg",
                    "--colDataFilePath", "colDataFilePathArg",
                    "--scalingFactorsPath", "estimatedScalingFactorsArg",
                    "--localCacheDir", "localCacheDirArg",
            ]

    cli.parse( inCommandLineArgsListMean )
    
    logger.info(f'command line arguments : {cli.args}')
    testLogCLI(logger, cli)

    # test useMedian 
    logger.info(f'########## test useMedian ')
    inCommandLineArgsListMedian = [
                    "--groupByGeneCountFilePath", "groupByGeneCountFilePathArg",
                    "--geneSignatureProfilesDataRootDir", "geneSignatureProfilesDataRootDirArg",
                    "--colDataFilePath", "colDataFilePathArg",
                    "--scalingFactorsPath", "estimatedScalingFactorsArg",
                    "--localCacheDir", "localCacheDirArg",
                    "--useMedian"
            ]
    cliMedian = CibersortSignatureMatrixFactoryCLI( 
                                             version=__version__ , 
                                             author=__author__ ,
                                             date=__date__, 
                                             update=__updated__)

    cliMedian.parse( inCommandLineArgsListMedian )

    logger.info(f'command line arguments : {cliMedian.args}')
    testLogCLI(logger, cliMedian)
 
    logger.info("END")

################################################################################
def main(inCommandLineArgsList=None):
    '''
    cd /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python
    (extraCellularRNA) aedavids@mustard $ python -m pipeline.dataFactory.cibersortMixtureMatrixFactory
    '''
    # we only configure logging in main module
    # loglevel = p.getProperty("LOG_LEVEL")
    loglevel = "INFO"
    # logFMT = p.getProperty("LOG_FMT")
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    logger = logging.getLogger(os.path.basename(__file__))

    # testCLI(logger)
    # logger.error('ERROR AEDWIP comment out testCLI!!!!!')
    # sys.exit(-1)

    cli = CibersortSignatureMatrixFactoryCLI( 
                                             version=__version__ , 
                                             author=__author__ ,
                                             date=__date__, 
                                             update=__updated__)

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.info(f'command line arguments : {cli.args}')
 
    geneSignatureProfilesDataRootDir = cli.args.geneSignatureProfilesDataRootDir
    groupByGeneCountFilePath = cli.args.groupByGeneCountFilePath
    colDataFilePath = cli.args.colDataFilePath
    scalingFactorsPath = cli.args.scalingFactorsPath
    localCacheDir = cli.args.localCacheDir
    outDir = cli.args.outDir
    if cli.args.useMedian:
        useMedian = True
    else :
        useMedian = False
    logger.info(f' useMedian : {useMedian}')


    csmf = CibersortSignatureMatrixFactory(
            geneSignatureProfilesDataRootDir=geneSignatureProfilesDataRootDir,
            groupByGeneCountFilePath=groupByGeneCountFilePath,
            colDataFilePath=colDataFilePath,
            estimatedScalingFactorsFilePath=scalingFactorsPath,
            localCacheDir=localCacheDir,
            outdir = outDir,
            testSize = None, 
            verbose = False,
            useMedian = useMedian
   )    

################################################################################
if __name__ == '__main__':
    main()
