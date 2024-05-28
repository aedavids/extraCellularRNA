#
# testOptimalSelectiveEnrichSignatureGeneConfig.py
# end to end test we can use to debug slurm scripts
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: testEnrichSignatureGeneConfig.py
#

import os
# AEDWIP_ENV = os.environ['AEDWIP_ENV']
# print(f'AEDWIP AEDWIP_ENV : {AEDWIP_ENV}')

ORIG_PYTHONPATH = os.environ['PYTHONPATH']
print(f'PYTHONPATH : {ORIG_PYTHONPATH}')
print(f'FILE: {__file__}')
print(f'PWD: {os.getcwd()}')

import ast

# https://docs.python.org/3/howto/logging.html#changing-the-format-of-displayed-messages
# https://stackoverflow.com/a/48996222
import logging

# import numpy as np
from numpy import nan
import pandas as pd
import pathlib as pl

import pprint as pp
import shutil
import unittest

from analysis.bestSignatureGeneConfig import BestSignatureGeneConfig
from analysis.optimalSelectiveEnrichSignatureGeneConfig import OptimalSelectiveEnrichSignatureGeneConfig
from analysis.utilities import findAllGenes, findFile, loadDictionary, saveSet
from pipeline.dataFactory.driver import _countExtraHeaderLines
from pipeline.dataFactory.driver import runSelectGenesOfInterest
from pipeline.dataFactory.upsetPlotDataFactory import UpsetPlotDataFactory
from plots.upsetPlots import UpsetPlot

################################################################################
class TestOptimalSelectiveEnrichSignatureGeneConfig(unittest.TestCase):
    '''
    TODO
    '''
    logger = logging.getLogger( os.path.basename(__file__) )

    # https://stackoverflow.com/a/14493895/4586180
    #maxDiff = None

    # python -m unittest discover .
    # we need to be able to find test file relative to module location
    # os.getcwd() is often the top of our source code tree
    # https://stackoverflow.com/a/57836145/4586180
    relativeRootPath = pl.Path(os.path.dirname(os.path.realpath(__file__)))
  
    # localCacheDir= "./data/tmp"
    localCacheDir= relativeRootPath.joinpath("data/testOptimalSelectiveEnrichSignatureGeneConfig/tmp")

    # The files produced by runSelectGenesOfInterest() should be found here
    filteredDESeqResultsDir = localCacheDir.joinpath('testDataSet-design-tilda-gender-category-padj-001-lfc-20-n-3')

    # expectedCandidateGenesDir = relativeRootPath.joinpath("data/testOptimalSelectiveEnrichSignatureGeneConfig/1vsAllExpected")
    # expectedCandidateGenes = [
    #     str(expectedCandidateGenesDir.joinpath("UVM_vs_all.results")),
    #     str(expectedCandidateGenesDir.joinpath("Vagina_vs_all.results")),
    #     str(expectedCandidateGenesDir.joinpath("Whole_Blood_vs_all.results"))
    # ]

    bestTopN = 3
    lfcThreshold = 2.0
    padjThreshold = 0.01
    
    GTEx_TCGADir = relativeRootPath.joinpath("data/testOptimalSelectiveEnrichSignatureGeneConfig/GTEx_TCGA1vsAll")
    historyGeneSetPath = localCacheDir.joinpath('historyGeneSetPath.txt')

    optimalSelectiveOut = localCacheDir.joinpath('SimulatedGTEx_TCGA-design-test_design-padj-001-lfc-20-n-5')

    upstreamOutDir          = localCacheDir.joinpath('SimulatedGTEx_TCGA-design-test_design-padj-001-lfc-20-n-3')
    upstreamUpsetPlotOutDir = localCacheDir.joinpath('upsetPlot.out')
    

    ################################################################################
    def setUp(self):
        '''
        runs before every test. Use to put system into a know state before each test
        '''
        self.logger.info("BEGIN")

        if os.path.exists(self.localCacheDir):
            self.logger.info(f"removing old localCacheDir  : {self.localCacheDir}")
            shutil.rmtree(self.localCacheDir)
            os.mkdir(self.localCacheDir)
        
        # signatureMatrixOutputFile = self.relativeRootPath.joinpath("data/geneSignatureProfiles/best/ciberSortInput/signatureGenes.tsv")
        # if os.path.exists(signatureMatrixOutputFile):
        #     self.logger.info(f"removing old output file : {signatureMatrixOutputFile}")
        #     os.remove(signatureMatrixOutputFile)           

        self.logger.info("END\n")

    ################################################################################
    @unittest.skip("skip template. it is not a real test")
    def testTemplate(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n")            

    ################################################################################
    def testOptimalSelectiveEnrich(self):
        self.logger.info("BEGIN")

        self._runUpstream()

        self._createHistoryGeneSet()

        sgc = self._getOptimalSelectiveEnrichGeneConfig()

        selectedGeneSetsDict, outFileList = self._runTestOptimalSelectiveEnrich(sgc)

        # for key,df in selectedGeneSetsDict.items():
        #     if key == 'Vagina_vs_all.results':
        #       there are two rows with same index values. to_dict() will drop one
        #         df = df.reset_index()
        #     print(f'\n\n{key} = {pp.pformat(df.to_dict())}\n')

        cols = ['name', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']

        expectedUVMDict = self.getExpectedResults(key=0)
        expecteUVMDF = pd.DataFrame(expectedUVMDict).loc[:, cols]

        expectedVaginaDict = self.getExpectedResults(key=1)
        expectedVaginaDF = pd.DataFrame(expectedVaginaDict).loc[:, cols]

        expectedWhole_BloodDict = self.getExpectedResults(key=2)
        expectedWhole_BloodDF = pd.DataFrame(expectedWhole_BloodDict).loc[:, cols]

        uvmDF = selectedGeneSetsDict['UVM_vs_all.results'].loc[:, cols]

        vaginaDF = selectedGeneSetsDict['Vagina_vs_all.results'].loc[:, cols]
        vaginaDF = vaginaDF.reset_index()
        selectCols = ~vaginaDF.columns.isin(['index'])
        vaginaDF = vaginaDF.loc[:, selectCols]

        whole_BloodDF = selectedGeneSetsDict['Whole_Blood_vs_all.results'].loc[:, cols]
        
        pd.testing.assert_frame_equal(expecteUVMDF, uvmDF)
        pd.testing.assert_frame_equal(expectedVaginaDF, vaginaDF)
        pd.testing.assert_frame_equal(expectedWhole_BloodDF, whole_BloodDF)

        # test
        # check degree1 dicts
        # check the results files
            
        self.logger.info("END\n")            

    ################################################################################
    def _createHistoryGeneSet(self) :
        '''
        includes all the genes that we have considered so far

        assume bestTopN was the pipeline head stage. Additional stages
        may have removed genes or add genes. Lets assume that stages that
        added genes provided a way to keep track of genes they considered

        lets not worry about this in our test
        '''
        self.logger.info("BEGIN")

        # It is easy to find all genes considered by bestTopN
        # they are all in the stages ouput directory
        pattern="*vs_all.results"
        upstreamResultsFiles = findFile( self.upstreamOutDir, pattern)

        historyGeneSet = set()
        for fPath in upstreamResultsFiles:
            numRowsToSkip = _countExtraHeaderLines(fPath)
            self.logger.info(f"fPath: {fPath}")
            deseqDF = pd.read_csv(fPath, skiprows=numRowsToSkip)
            genes = deseqDF.loc[:, "name"].to_list()
            historyGeneSet = historyGeneSet.union(genes)

        saveSet(self.historyGeneSetPath, historyGeneSet)

        self.logger.info("END")
    
    ################################################################################
    def _getOptimalSelectiveEnrichGeneConfig(self) -> OptimalSelectiveEnrichSignatureGeneConfig:
        '''
        todo
        '''

         # Vagina has 1 degree1 gene, UVM does not have any
        categories = [ 'Vagina', 'UVM']

        # all upstream have 3 degree 1 genes
        maxNumberOfGenes = 4

        # this is the directory our pipeline head task, bestTopN started with
        deseqResultsDir = self.GTEx_TCGADir

        upStreamIntersectionDictionaryPath = self.upstreamUpsetPlotOutDir.joinpath('test_title.intersection.dict')
        sgc = OptimalSelectiveEnrichSignatureGeneConfig(
                dataSetName="SimulatedGTEx_TCGA", 
                design="test_design", 
                padjThreshold=self.padjThreshold, 
                lfcThreshold=self.lfcThreshold, 
                windowLength=4, # our test data only has 9 rows range will be 3 ... 3 + 4
                localCacheRootPath=str(self.localCacheDir), 
                title="test_title",
                deseqResultsDir=deseqResultsDir,
                categories=categories,
                startIdx=self.bestTopN,
                maxNumberOfGenes=maxNumberOfGenes,
                upStreamIntersectionDictionaryPath=upStreamIntersectionDictionaryPath,
                historicGeneSetPath =self.historyGeneSetPath,
                )
        
        return sgc

    ################################################################################
    def _runTestOptimalSelectiveEnrich(self,  
                                       sgc : OptimalSelectiveEnrichSignatureGeneConfig)-> tuple[dict, list]:
        '''
        TODO
        '''
        self.logger.info('BEGIN')

        pattern="*vs_all.results"
        upstreamResultsFile = findFile(self.upstreamOutDir, pattern)
        self.logger.info(f'pattern: \n{pattern}')

        selectedGeneSetsDict, outFileList = runSelectGenesOfInterest(
                                                signatureGeneConfig=sgc, 
                                                candidateSignatureFileList=upstreamResultsFile ,
                                            ) 
        
        self.logger.info('END')
        return (selectedGeneSetsDict, outFileList)

    ################################################################################
    def _runUpstream(self) :
        '''
        this is our pipeline's head stage
        run bestN 
        '''
        self.logger.info('BEGIN')
        sgc = BestSignatureGeneConfig(
                                    dataSetName="SimulatedGTEx_TCGA",
                                    design="test_design",
                                    padjThreshold=self.padjThreshold,
                                    lfcThreshold=self.lfcThreshold,
                                    n=self.bestTopN,
                                    localCacheRootPath=str(self.localCacheDir),
                                    title="test_title"
                                    )         
        
        pattern="*vs_all.results"
        GTEx_TCGA1vsAllResultsFile = findFile(self.GTEx_TCGADir, pattern)
        self.logger.info(f'GTEx_TCGA1vsAllResultsFile: \n{GTEx_TCGA1vsAllResultsFile}')
        selectedGeneSetsDict, outFileList = runSelectGenesOfInterest(
                                            signatureGeneConfig=sgc, 
                                            candidateSignatureFileList=GTEx_TCGA1vsAllResultsFile ,
                                            ) 
        
        upsetFactory = UpsetPlotDataFactory()
        upsetPlotDataDF, geneSetsDict = upsetFactory.createUpsetPlotDataFromSetDict(selectedGeneSetsDict)
        self.logger.info(f'upsetPlotDataDF:\n{pp.pformat(upsetPlotDataDF)}')

        upsetPlot = UpsetPlot(signatureGeneConfig=sgc, 
                        geneSetsDict=geneSetsDict, 
                        geneSetsUpsetPlotData=upsetPlotDataDF)
        
        intersectionDict = upsetPlot.findIntersectionElements()  
        upsetPlot.saveInteresection(self.upstreamUpsetPlotOutDir, intersectionDict)

        self.logger.info('END')

    ################################################################################   
    def getExpectedResults(self, key):
        '''
        todo
        UVM = 0
        Vagina = 1
        Whole_Blood = 2
        '''
        expected = [
                {'baseMean': {0: 26.107766460079,
                                                            1: 18.4753729001426,
                                                            2: 12.8920299186668,
                                                            4: 11.8549907875931},
                                                'lfcSE': {0: 0.953575988058649,
                                                        1: 0.849704091672415,
                                                        2: 0.841781446504066,
                                                        4: nan},
                                                'log2FoldChange': {0: -24.5945461542676,
                                                                    1: -20.2185391172986,
                                                                    2: -23.7393571218277,
                                                                    4: -24.9029360519086},
                                                'name': {0: 'UVM_V', 1: 'AC114812.1', 2: 'UVM_V_W.1', 4: 'UVM_AAA'},
                                                'padj': {0: 1.08832945564078e-143,
                                                        1: 1.7198084154408898e-122,
                                                        2: 1.15291027996521e-171,
                                                        4: 1.0730900557478e-144},
                                                'pvalue': {0: 1.09281887371659e-146,
                                                            1: 3.7805167621231e-125,
                                                            2: 5.63188912271324e-175,
                                                            4: nan},
                                                'stat': {0: -25.7919100965816,
                                                        1: -23.7948002315769,
                                                        2: -28.2013309041411,
                                                        4: nan}},
                
                # {'baseMean': {0: 2958.09826948922, 1: 510.344123086478, 2: 1323.93436332429},
                #                         'lfcSE': {0: 0.380337343444453, 1: nan, 2: 0.421037362801773},
                #                         'log2FoldChange': {0: -7.63510684000441,
                #                                             1: -7.93491125619565,
                #                                             2: -7.62072609459715},
                #                         'name': {0: 'C1orf61', 1: 'OLIG2', 2: 'OLIG1'},
                #                         'padj': {0: 3.17969371668197e-86,
                #                                 1: 7.01524601037036e-50,
                #                                 2: 2.355818959958549e-70},
                #                         'pvalue': {0: 1.23155292517976e-89, 1: nan, 2: 3.19357962258532e-73},
                #                         'stat': {0: -20.0745653078883, 1: nan, 2: -18.0998808369057}},

                {'baseMean': {0: 2958.09826948922,
                            1: 1709.82218160823,
                            2: 1323.93436332429,
                            3: 510.344123086478},
                #'index': {0: 0, 1: 1, 2: 2, 3: 1},
                'lfcSE': {0: 0.380337343444453,
                        1: 0.39938008608752,
                        2: 0.421037362801773,
                        3: nan},
                'log2FoldChange': {0: -7.63510684000441,
                                    1: -8.02201914088103,
                                    2: -7.62072609459715,
                                    3: -7.93491125619565},
                'name': {0: 'C1orf61', 1: 'ARPP21', 2: 'OLIG1', 3: 'OLIG2'},
                'padj': {0: 3.17969371668197e-86,
                        1: 2.71054099468817e-86,
                        2: 2.355818959958549e-70,
                        3: 7.01524601037036e-50},
                'pvalue': {0: 1.23155292517976e-89,
                            1: 9.74852900208769e-90,
                            2: 3.19357962258532e-73,
                            3: nan},
                'stat': {0: -20.0745653078883,
                        1: -20.0861771037905,
                        2: -18.0998808369057,
                        3: nan}},
                                       
                {'baseMean': {0: 985081.569770207, 1: 939295.036780794, 2: 384428.792719921},
                                        'lfcSE': {0: 0.11089940519715, 1: 0.106055647081676, 2: 0.107194068493031},
                                        'log2FoldChange': {0: 11.3779865954891,
                                                            1: 11.3962888262948,
                                                            2: 11.3333149379708},
                                        'name': {0: 'V_W', 1: 'W_CCC', 2: 'UVM_V_W'},
                                        'padj': {0: 0.0, 1: 0.0, 2: 0.0},
                                        'pvalue': {0: 0.0, 1: 0.0, 2: 0.0},
                                        'stat': {0: 102.597363577037, 1: 107.455747429632, 2: 105.727071444328}}
        ]

        return expected[key]

# ################################################################################   
# def _get1vsAllUpstreamResults(relativeRootPath : pl.Path) -> list[str]:
#         '''
#         upstream task is testEnrichSignatureGeneConfig
#         these are example of DESEq results file produced by 1vsAll.
#         they have multiline headers and thousands of rows (they have not been filtered)
#         '''
#         retList = [
#             str(relativeRootPath.joinpath("data/testEnrichSignatureGeneConfig/1vsAll/UVM_vs_all.results")),
#             str(relativeRootPath.joinpath("data/testEnrichSignatureGeneConfig/1vsAll/Vagina_vs_all.results")),
#             str(relativeRootPath.joinpath("data/testEnrichSignatureGeneConfig/1vsAll/Whole_Blood_vs_all.results")),
#         ]
#         return retList

   ################################################################################
    @unittest.skip("skip. it is not a real test")
    def testReverseEngineerBest500FindAllDegree1_wl500_sh(self):
        '''
        Best500FindAllDegree1_wl500 what does this do?

        seems to just return best500
        '''
        self.logger.info("BEGIN")

        # values from /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/
        # best500FindAllDegree1_wl500/training/best500FindAllDegree1_wl500.sh.log
        colDataPath = '/private/groups/kimlab/GTEx_TCGA/groupbyGeneTrainingSets/GTEx_TCGA_TrainColData.csv'
        countDataPath = '/private/groups/kimlab/GTEx_TCGA/groupbyGeneTrainingSets/GTEx_TCGA_TrainGroupby.csv' 
        deseqResultsDir = '/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500GTEx_TCGA/training/best500GTEx_TCGA.sh.out/GTEx_TCGA-design-tilda_gender_category-padj-0001-lfc-20-n-500'
        estimatedScalingFactors = '/private/groups/kimlab/GTEx_TCGA/1vsAll/estimatedSizeFactors.csv' 
        findModule = 'analysis.createOptimalSelectiveEnrichSignatureGeneConfig' 
        #outDir = '/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500FindAllDegree1_wl500/training/./best500FindAllDegree1_wl500.sh.out'
        #--vargs 
        design = 'tilda_gender_category '
        padjThreshold = 0.001 
        lfcThreshold = 2.0 
        dataSetName = 'GTEx_TCGA '
        windowLength = 500 
        title = 'best500_findAllDegree1_wl500'
        # --localCacheRoot /private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500FindAllDegree1_wl500/training/./best500FindAllDegree1_wl500.sh.out 
        classes = 'Brain_Amygdala LIHC TGCT COAD Lung Esophagus_Gastroesophageal_Junction Brain_Substantia_nigra Kidney_Cortex THCA Adrenal_Gland Breast_Mammary_Tissue ACC KIRP LGG Uterus Artery_Coronary BRCA HNSC Nerve_Tibial Brain_Putamen_basal_ganglia KICH Brain_Caudate_basal_ganglia Brain_Nucleus_accumbens_basal_ganglia Adipose_Subcutaneous Brain_Hippocampus CHOL Cells_Cultured_fibroblasts Cells_EBV-transformed_lymphocytes PCPG Skin_Not_Sun_Exposed_Suprapubic Adipose_Visceral_Omentum Colon_Sigmoid Brain_Cerebellum SKCM Brain_Cerebellar_Hemisphere DLBC Pituitary Small_Intestine_Terminal_Ileum Vagina Brain_Frontal_Cortex_BA9 Ovary MESO Heart_Atrial_Appendage STAD Colon_Transverse Esophagus_Muscularis Whole_Blood Artery_Aorta Pancreas THYM Cervix_Endocervix BLCA Brain_Anterior_cingulate_cortex_BA24 Brain_Hypothalamus Skin_Sun_Exposed_Lower_leg LUSC LUAD PAAD CESC UCS GBM Brain_Cortex Heart_Left_Ventricle Brain_Spinal_cord_cervical_c-1 READ Minor_Salivary_Gland Stomach PRAD Testis OV Esophagus_Mucosa Prostate Artery_Tibial Muscle_Skeletal UCEC Bladder Spleen KIRC Thyroid SARC Liver ESCA UVM'
        categories = classes
        localCacheRootPath="./"
        historicGeneSetPath = '/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500FindAllDegree1_wl500/training/./best500FindAllDegree1_wl500.sh.out/historicGeneSet.txt'
        maxNumberOfGenes = 10000 
        resultsDir = '/private/groups/kimlab/GTEx_TCGA/1vsAll'
        startIdx = 500 
        upStreamIntersectionDictionaryPath = '/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category/best500GTEx_TCGA/training/best500GTEx_TCGA.sh.out/upsetPlot.out/best500.intersection.dict'
            
        osesgc = OptimalSelectiveEnrichSignatureGeneConfig(
                                dataSetName, 
                                design, 
                                padjThreshold, 
                                lfcThreshold, 
                                windowLength, 
                                localCacheRootPath, 
                                title,
                                deseqResultsDir,
                                #upstreamDeseqResultsDir : str, do we need this? _step1 will pass this to findGenes()
                                categories,
                                startIdx,
                                maxNumberOfGenes,
                                upStreamIntersectionDictionaryPath,
                                historicGeneSetPath
        )

        fileName = f'ESCA_vs_all.results'
        ESCAPath = f'{deseqResultsDir}/{fileName}'
        escaDESeqResultsDF = pd.read_csv(ESCAPath) #, skiprows=numRowsToSkip)
        osesgc.findGenes(escaDESeqResultsDF, fileName)

        self.logger.info("END\n")            

################################################################################
if __name__ == "__main__":
    # we only configure logging in main module
    # loglevel = p.getProperty("LOG_LEVEL")
    loglevel = "INFO"
    # logFMT = p.getProperty("LOG_FMT")
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
