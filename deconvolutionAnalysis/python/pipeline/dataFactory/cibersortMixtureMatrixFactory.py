#
# cibersortMixtureFactory.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/createCibersortMixtureMatrix.ipynb
#

# CibersortMixtureFactoryCommandLine display the doc string 
'''
    Converts a transcript aggregated by gene loci counts to ciber sort mixture matrix format. 
    See extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/test/data/ 
    for example input data file formats.
'''

import logging
import numpy as np
import os
import pandas as pd
import pathlib as pl
from pipeline.dataFactory.utilities import loadCache
from pipeline.dataFactory.cibersortMixtureMatrixFactoryCLI import CibersortMixtureFactoryCommandLine

__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-10-04'
__updated__ = '2023-10-04'

###############################################################################
class CibersortMixtureFactory( object ):
    '''
    ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/createCibersortMixtureMatrix.ipynb

    
    public functions. See doc string for details
        __init__()
        getExpectedFractionsDF()
        getLabledMixtureDF()  
        randomizeMixture(seed=None)
        saveMixtureAndExpectedFractions()
        saveRandomizedMixture()
    '''
    logger = logging.getLogger(__name__)

    ################################################################################    
    def __init__(self, 
                 signatueGeneFilePath : str,
                 groupByGeneCountFilePath : str,
                 colDataFilePath : str,
                 scalingFactorsPath : str,
                localCacheDir : str, 
                ) :
        '''
        arguments:
             arguments:
            signatueGeneFilePath : str
                path to a tsv file created by createCiberSortGeneSignatureMatrix.ipynb
                Each row in this file coresponds to a 1vsAll DESeq result for a specific gene
                
            groupByGeneCountFilePath : str
                path to a csv file with gene counts. 
                ex. 'groupbyGeneTrainingSets/GTEx_TCGA_TrainGroupby.csv'
                
            colDataFilePath : str
                path to a csv file containing sample meta data in DESeq format
            
            scalingFactorsPath : str
                path to a csv file of DESeq estimated scaling factors used 
                to adjust each sample to account for libaray size and library composition

            localCacheDir : str
                path to local file cache
                example: /scratch/aedavids/tmp

                when reading files will first check localCache. If not in cache will copy
                ensures file are localized.
                         
        '''
        self.localCacheDir = localCacheDir
        self.logger.info("\n localCacheDir:\n{}".format(self.localCacheDir))

        #self.signatueGeneFilePath = loadCache( signatueGeneFilePath, self.localCacheDir)
        self.signatueGeneFilePath =  signatueGeneFilePath
        self.logger.info("\n gene signature file:\n{}".format(self.signatueGeneFilePath))

        self.groupByGeneCountFilePath = loadCache(groupByGeneCountFilePath, self.localCacheDir)
        self.logger.info("\n groupByGeneCountFilePath:\n{}".format(self.groupByGeneCountFilePath))
        
        self.colDataFilePath = loadCache( colDataFilePath, self.localCacheDir)
        self.logger.info("\n colDataFilePath:\n{}".format(self.colDataFilePath))

        self.scalingFactorsPath = loadCache( scalingFactorsPath, self.localCacheDir)
        self.logger.info("\n scalingFactorsPath:\n{}".format(self.scalingFactorsPath))
     
        #self.outdir = outdir
        
        self.metaColList = ['sample_id', 'participant_id', 'category', 'gender', 'age', 'dataSet']

        self.geneSignatureDF = None
        self.groupedByGeneDF = None
        self.colDataDF = None
        self.scalingFactorDF = None
        self._load()

        self._normalize()
        
        self.geneList = None
        self.mixtureDF = None
        self._select()
                
        self.labeledMixtureDF = None
        self._createLabeledMixtue()
        
        self.expectedFractionsDF = None 
        self._createFractions()
    
    ################################################################################
    def getExpectedFractionsDF(self):
        return self.expectedFractionsDF
    
    ################################################################################
    def getLabeledMixtureDF(self):
        return self.labeledMixtureDF
    

    ################################################################################
    def randomizeMixture(self, seed=None):
        '''
        returns a randomly shuffled copy of the mixtureDF. Randomizing data remove all 
        information, creating a good worst case base line you can use to evaluate
        models agains

        arguments:
            seed:
                integer
                set if want the psudo random generator to return the same sequence of
                results. Use for testing purpose purposes only


        ref: 
            - https://github.com/aedavids/lab3RotationProject/blob/master/src/test/testDataFactory.py)]
            - https://github.com/aedavids/lab3RotationProject/blob/master/src/DEMETER2/dataFactory.py
        '''
        if seed:
            np.random.seed(seed)
        else:
            epochTime = int(time.time())
            np.random.seed(epochTime)


        # make a copy to ensure no side effect
        # remove the gene names and meta data. we only want to 
        # shuffle the count values
        copyDF = self.labeledMixtureDF.copy()   
        
        cols = copyDF.columns.to_list()
    
        removedMetaCols = []
        for colName in self.metaColList:
            if colName in cols:
                cols.remove(colName)
                removedMetaCols.append(colName)

    
        copyDF = copyDF.loc[:, cols]        
        removedMetaDF =  self.labeledMixtureDF.loc[:, removedMetaCols] 

        numRows, numCols = copyDF.shape   
        valuesNP = copyDF.values
        for r in range(numRows):
            randomRowIdx = np.random.permutation(numCols)
            #self.logger.info("randomRowIdx:{}".format(randomRowIdx))
            # use numpy fancy indexing
            valuesNP[r,:] = valuesNP[r, randomRowIdx]

        # create data frame from the scrambled valuesNP
        randomDF = pd.DataFrame(valuesNP, columns=copyDF.columns )

        # add the meta data back
        byColumns=1 # column bind
        retDF = pd.concat([randomDF, removedMetaDF], axis=byColumns)
    
        return retDF
    
    
    ################################################################################
    def saveMixtureAndExpectedFractions(self, outDir, prefixStr=None) -> tuple[pl.Path, pl.Path]:
        '''
        saves in a cibersort's expected format

        returns : tuple[pl.Path, pl.Path]
                (mixturePath, expectedFractionsPath)
        '''
#         base) $ cat mixture.txt 
#         sampleTitle	S1	S2	S3	S4	S5	S6
#         G1	1.0	0.0	0.0	1.0	1.0	0.0
#         G2	1.0	0.0	0.0	1.0	1.0	0.0
        ciberSortFmtDF = self._convertToMixtureToCiberSortFmt(self.labeledMixtureDF)
        mixturePath = self._save( outDir, ciberSortFmtDF, fileName="mixture.txt", prefixStr=prefixStr)
        
        expectedFractionsPath = self._save( outDir, self.expectedFractionsDF, fileName="expectedFractions.txt", prefixStr=prefixStr)

        return (mixturePath, expectedFractionsPath)
    
    
    ################################################################################
    def _convertToMixtureToCiberSortFmt(self, df):
        '''
        TODO
        '''
        metaList = self.metaColList.copy()
        metaList.remove('sample_id')
        dataCols = df.columns.to_list()
        for m in metaList:
            dataCols.remove(m)

        #self.logger.info(dataCols)
        mixtureDF = df.loc[:, dataCols].copy()
        #display(mixtureDF.shape)
        #mixtureDF.head()

        #self.logger.info("\n rename and set index")
        mixtureDF = mixtureDF.rename(columns={"sample_id":"sampleTitle"})
        mixtureDF = mixtureDF.set_index('sampleTitle')
        #self.logger.info(mixtureDF.shape)
        #display( mixtureDF.head() )

        #self.logger.info("\n transpose")
        mixtureDF = mixtureDF.transpose()
        #display( mixtureDF.head() )
        
        #self.logger.info("\n reset index")
        mixtureDF = mixtureDF.reset_index()
        # the original index data is now a column with name 'index'
        mixtureDF = mixtureDF.rename(columns={"index":"sampleTitle"})
        #display( mixtureDF.head() ) 
        
        return mixtureDF

    
    ################################################################################
    def saveRandomizedMixture(self, outDir, randomizedDF, prefixStr=None) -> pl.Path:
        '''
        randomization removes all information from the mixture. It can be used to 
        establish a base line to compared trained models againt

        arguments:
            randomizedDF = saveRandomizedMixture()

        saves in a cibersort's expected format

        returns : pl.Path
        '''
        ciberSortFmtDF = self._convertToMixtureToCiberSortFmt(randomizedDF)        
        retPath = self._save( outDir, ciberSortFmtDF, fileName="randomizedMixture.txt", prefixStr=prefixStr) 

        return retPath  
    
    ################################################################################
    def _save(self, outDir, df, fileName, prefixStr=None) -> pl.Path:
        '''
        saves in a cibersort's expected format
        
        arguments:
            fileName
                string. example 'mixture.txt'
        '''
        dataOutdir = pl.Path(outDir)

        if prefixStr:
            path = dataOutdir.joinpath(prefixStr + "_" + fileName)
        else:
            path = dataOutdir.joinpath(fileName)   

        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False, sep="\t")
        self.logger.info("\n saved to: {}".format(path))    

        return path      
        
    ################################################################################
    def _load(self):
        self.geneSignatureDF = pd.read_csv(self.signatueGeneFilePath, sep="\t" )
        self.logger.info("\ngeneSignatureDF.shape:{}".format(self.geneSignatureDF.shape))
        self.logger.debug("geneSignatureDF.iloc[0:3, :]")
        self.logger.debug(self.geneSignatureDF.iloc[0:3, :])
        
        self.groupedByGeneDF = pd.read_csv(self.groupByGeneCountFilePath, sep=",")
        self.logger.info("\ngroupByGeneDF.shape:{}".format(self.groupedByGeneDF.shape))
        self.logger.debug("groupByGeneDF.iloc[0:3, 0:3]")
        self.logger.debug(self.groupedByGeneDF.iloc[0:3, 0:3])
        
        self.colDataDF = pd.read_csv(self.colDataFilePath, sep=",")
        self.logger.debug("colDataDF.iloc[0:3, :]")
        self.logger.debug(self.colDataDF.iloc[0:3, :])
        
        self.scalingFactorDF = pd.read_csv(self.scalingFactorsPath, sep=",")
        self.logger.debug("scalingFactorDF.iloc[0:3, :]")
        self.logger.debug(self.scalingFactorDF.iloc[0:3, :])
        
    ################################################################################
    def _normalize(self):
        #aedwip we can not multiple has geneIdcol
        cols = list(self.groupedByGeneDF.columns)
        cols.remove('geneId')
        self.logger.info(f'groupedByGeneDF \n{self.groupedByGeneDF}')

        self.logger.info(f'groupedByGeneDF.shape : {self.groupedByGeneDF.shape}')
        self.logger.info(f'scalingFactorDF.shape : {self.scalingFactorDF.shape}')
        # aedwip = self.groupedByGeneDF * self.scalingFactorDF.transpose().values
        # self.logger.info(f'aedwip \n{aedwip}')

        countsDF = self.groupedByGeneDF.loc[:, cols]
        normalizedCountsDF = countsDF * self.scalingFactorDF.transpose().values
        self.groupedByGeneDF.loc[:, cols] = normalizedCountsDF
        self.logger.debug(f'normalized \n{self.groupedByGeneDF}')

    ################################################################################
    def _select(self):
        self.geneList = self.geneSignatureDF.loc[:, "name"].to_list()
        #oneVsAllDF = oneVsAllDF[ oneVsAllDF.loc[:,"name"].isin( self.geneListsorted) ]
        #oneVsAllDF = oneVsAllDF.sort_values( by=["name"] )
        
        gbDF = self.groupedByGeneDF
        df = gbDF[ gbDF.loc[:,"geneId"].isin(self.geneList) ]
        # sort makes debug easier
        df = df.sort_values( by=["geneId"] )

        # do not change sort order. It shold already match colData order
        
        # rename to match cibersort expected format
        df = df.rename(columns={'geneId':"name"})
        
        # set index will make join after traspose easier
        # if we do not set the index, the 'name' column will become a row
        # instead of the column names
        df = df.set_index('name')
        
        self.mixtureDF = df
        self.logger.debug("\n mixtureDF.iloc[0:3, 0:3]")
        self.logger.debug(self.mixtureDF.iloc[0:3, 0:3])
        

    ################################################################################
    def _createLabeledMixtue(self):
        '''
        transpose,
        scale,
        join mixtureDF with colData
        '''
        # # normalize counts      
        self.logger.info(f'mixtureDF.shape : {self.mixtureDF.shape}')
        self.logger.info(f'mixtureDF :\n{self.mixtureDF}')
       

        transposeGroupByDF = self.mixtureDF.transpose(copy=True)
        self.logger.info(f'transposeGroupByDF.shape : {transposeGroupByDF.shape}')
        self.logger.info(f'transposeGroupByDF :\n{transposeGroupByDF}')

        self.logger.info(f'colDataDF.shape : {self.colDataDF.shape}')
        self.logger.info(f'colDataDF :\n{self.colDataDF}')

        self.labeledMixtureDF =  pd.merge(left=transposeGroupByDF, 
                right=self.colDataDF, 
                how='inner', 
                left_index=True, #left_on="name", 
                right_on="sample_id")

        self.logger.info("\n labeledMixtureDF.iloc[0:3, 0:3]")
        self.logger.info(self.labeledMixtureDF.iloc[0:3, 0:3])  
        
        #metaColList = ['sample_id', 'participant_id', 'category', 'gender', 'age', 'dataSet']
        self.logger.debug("\n labeledMixtureDF.iloc[0:3, {}]".format(self.metaColList))
        self.logger.debug(self.labeledMixtureDF.loc[0:3, self.metaColList])        
        
    ################################################################################    
    def _createFractions(self):
        '''
        TODO
        '''
        
        df = self.labeledMixtureDF.loc[:,["category"]].drop_duplicates()\
                    .sort_values(by='category')
        listOfTypes = df["category"].values.tolist()
        self.logger.info("\n number of types: {}".format(len(listOfTypes)))
        
        # create an empty data frame, and set columns
        self.expectedFractionsDF = pd.DataFrame(columns=['sample_id'] + listOfTypes)
        numTypes = len(listOfTypes)

        for index, row in self.labeledMixtureDF.iterrows():
            sample_id = row['sample_id']
            category = row['category']
            #self.logger.info("sample_id:{} category:{}".format(sample_id, category))
            idx = listOfTypes.index( category)
            linearCombination = np.zeros(numTypes)
            linearCombination[idx] = 1.0
            # use comprehension to expand numpy array into format that works with pandas
            # this method of append seems slow, it okay for GTEX_TCGA
            # for large data sets consider
            # https://stackoverflow.com/a/48287388
            self.expectedFractionsDF.loc[ len(self.expectedFractionsDF.index) ] = [sample_id]\
                            + [i for i in linearCombination]

        # pre pend the meta data to the extected fractions data frame
        # would have to write a lot more buggy code if we did this in for loop
        #metaColList = ['sample_id', 'participant_id', 'category', 'gender', 'age', 'dataSet']
        metaDF = self.labeledMixtureDF.loc[:, self.metaColList]
        self.expectedFractionsDF = pd.merge(left=metaDF, 
                                            right=self.expectedFractionsDF, 
                                            how='inner', 
                                            left_on="sample_id",
                                            right_on="sample_id")

        self.logger.debug(f"\n expectedFractionsDF\n{self.expectedFractionsDF}")

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

    logger = logging.getLogger(__file__)
    cli = CibersortMixtureFactoryCommandLine( 
                                             version=__version__ , 
                                             author=__author__ ,
                                             date=__date__, 
                                             update=__updated__)

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.info(f'command line arguments : {cli.args}')
 
    signatueGeneFilePath = cli.args.signatueGeneFilePath
    groupByGeneCountFilePath = cli.args.groupByGeneCountFilePath
    colDataFilePath = cli.args.colDataFilePath
    scalingFactorsPath = cli.args.scalingFactorsPath
    localCacheDir = cli.args.localCacheDir
    outDir = cli.args.outDir
    prefix = cli.args.prefix

    csmf = CibersortMixtureFactory(
                    signatueGeneFilePath,
                    groupByGeneCountFilePath, 
                    colDataFilePath,
                    scalingFactorsPath,
                    localCacheDir
                )

    mixturePath, expectedFractionsPath = csmf.saveMixtureAndExpectedFractions(outDir=outDir, 
                                                                              prefixStr=prefix)
    logger.info(f'mixture saved to: {mixturePath}')
    logger.info(f'expectedFractions save to: {expectedFractionsPath}')

    randomizedDF = csmf.randomizeMixture(seed=42)
    randomizedPath = csmf.saveRandomizedMixture(outDir=outDir, 
                                                randomizedDF=randomizedDF,
                                                prefixStr=prefix)

################################################################################
if __name__ == '__main__':
    main()
