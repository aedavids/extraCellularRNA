#
# pipeline.py
# main(): runs deconvolution hyper parameter tunning pipeline
# see pipelineCLITest.sh
# ref: extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/test/testIntegration.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
# 11/13/23
#

# upstreamPiplineCLI CommandLine display the doc string 
'''
    prepares DESeq and transcript count files so that we run CIBERSORTx
        1. select genes of interest
        2. create upset plots and find interesection elements
        3. create CIBERSORTx signatureGenes matrix
        4. create CIBERSORTx mixture matrix    
    
    See extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/test/testIntegration.py

'''

import glob
import importlib
import logging
import os

from pipeline.dataFactory.cibersortMixtureMatrixFactory import CibersortMixtureFactory
from pipeline.dataFactory.cibersortSignatureMatrixFactory import CibersortSignatureMatrixFactory
from pipeline.dataFactory.driver import runSelectGenesOfInterest
from pipeline.dataFactory.signatureGeneConfig import SignatureGeneConfiguration
from pipeline.dataFactory.upsetPlotDataFactory import UpsetPlotDataFactory

from pipeline.upstreamPipelineCLI import UpstreamPipelineCommandLine
#use this import if in visual studio code from pipelineCLI import PipelineCommandLine

import pprint as pp

import sys
from plots.upsetPlots import UpsetPlot

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-11-13'
__updated__ = '2023-11-13'


################################################################################
def _step1(logger           : logging.Logger,
           deseqResultsDir  : str,
           findModule       : str,
           outDir           : str,
           vargs            : list
           ) -> tuple[SignatureGeneConfiguration, dict, list ]:
    '''
    step 1.
        - create user defined SignatureGeneConfiguration obj
        - run SignatureGeneConfiguration.findGenes() on all candidate results file
    
    returns 
        (sgc, selectedGeneSetsDict, outFileList)
    '''
    logger.info(f'BEGIN')
    try:
        logger.info('step 1: create user defined SignatureGeneConfiguration obj')
        myModule = importlib.import_module(findModule, package=None)
        sgc = myModule.createSignatureGeneConfig(outDir, vargs)
        logger.info(f'sgc : {sgc}')
        # AEDWIP logger statement only works if '--findModule  pipeline.dataFactory.test.exampleCreateSignatureGeneConfig'
        #retStr = sgc.testDynamicalyloadedModules()
        #logger.info(f"sgct.testDynamicalyloadedModules() : {retStr}")

    except ImportError as e:
        emsg = f'unable to dynamically import {findModule} exception: {e}'
        logger.error(emsg)
        sys.exit(1)  

    #
    # get a list of deseq results file we want to select gene signatures from
    #
    deseqResultFiles = os.path.join(deseqResultsDir, "*.results")
    deseqResultFileList = glob.glob(deseqResultFiles)
    logger.warning(f'deseqResultFileList : \n{pp.pformat(deseqResultFileList)}')


    #
    # runSelectGenesOfInterest
    # filters/selects the genes of interest from each of the 1vsAll DESeq  results file
    # the selectedGeneSetsDict contains a key/value for each file
    # the outFileList is a file with a list of path to filtered DESeq results file.
    # ie only DESeq results we are interested
    #

    (selectedGeneSetsDict, outFileList) = runSelectGenesOfInterest(
                                    signatureGeneConfig=sgc, 
                                    candidateSignatureFileList=deseqResultFileList 
                                    ) 

    logger.info("END")
    return (sgc, selectedGeneSetsDict, outFileList) 

################################################################################
def _step2(logger           : logging.Logger,
           outDir : str,
           selectedGeneSetsDict : dict,
           sgc : SignatureGeneConfiguration,
           ):
    '''
    # create upset plots and find interesection elements
    '''
    logger.info("BEGIN")
    logger.info('step 2 create plots and find intersection elements')
    upsetFactory = UpsetPlotDataFactory()
    upsetPlotDataDF, geneSetsDict = upsetFactory.createUpsetPlotDataFromSetDict(selectedGeneSetsDict)
    logger.info(f'upsetPlotDataDF:\n{pp.pformat(upsetPlotDataDF)}')

    upsetPlot = UpsetPlot(signatureGeneConfig=sgc, 
                    geneSetsDict=geneSetsDict, 
                    geneSetsUpsetPlotData=upsetPlotDataDF)
    
    #
    # create a plot of all intersections that are not empty
    # 
    fig, pltDict = upsetPlot.configurablePlot(show_counts=True)
    upsetPlotOutDir = os.path.join(outDir, "upsetPlot.out")
    plotPath = upsetPlot.savePlot( outdir=upsetPlotOutDir,
                                    extraFileNameParts = 'allDegrees')
    logger.warning(f'plotPath : {plotPath}')

    intersectionDict = upsetPlot.findIntersectionElements()  
    intersectionDictPath = upsetPlot.saveInteresection(upsetPlotOutDir, intersectionDict) 
    logger.warning(f'intersectionDictPath : {intersectionDictPath}')

    #
    # the plot of all intersections can be difficult to use
    # for each degree create a separate plot
    # the degree of an intersection is the number of sets it was composed of
    # the degree for A.intersections(B) is 2. 
    #
    degrees = upsetPlot.findDegrees(intersectionDict)
    logger.info(f'degrees: {degrees}')

    for degree in degrees:
        min_degree = degree
        max_degree = degree
        fig, pltDict = upsetPlot.configurablePlot(show_counts=True, min_degree=min_degree, max_degree=max_degree)
        vargs = f"min_degree={min_degree},_max_degree={max_degree}"
        plotPath = upsetPlot.savePlot( outdir=upsetPlotOutDir,
                                extraFileNameParts=vargs)
        logger.info(f'plotPath : {plotPath}')    

    logger.info("END")

################################################################################
def main(inCommandLineArgsList=None):
    '''
    ref: pipeline/dataFactory/test/testInjection.py 

    ```
    # run using unit test data
    cd /private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python
    $ see pipelineCLITest.sh
    ```
    '''
    # we only configure logging in main module
    #loglevel = "WARN"
    loglevel = "INFO"
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    logger = logging.getLogger("__name__")
    logger.warning("BEGIN")

    #
    # always log run time env to make debug easier
    #
    ORIG_PYTHONPATH = os.environ['PYTHONPATH']
    logger.warning(f'PYTHONPATH : {ORIG_PYTHONPATH}')
    logger.warning(f'FILE: {__file__}')
    logger.warning(f'PWD: {os.getcwd()}')

    #
    # parse command line arguments
    #
    cli = UpstreamPipelineCommandLine( 
                            version=__version__ , 
                            author=__author__ ,
                            date=__date__, 
                            update=__updated__)

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    colDataPath     = cli.args.colDataPath
    countDataPath    = cli.args.countDataPath
    deseqResultsDir = cli.args.deseqResultsDir
    esfPath         = cli.args.estimatedScalingFactors
    findModule      = cli.args.findModule
    outDir          = cli.args.outDir
    vargs           = cli.args.vargs

    logger.info(f"vargs: {vargs}")

    #
    # get upset plot input files
    #
    sgc, selectedGeneSetsDict, outFileList = _step1(logger, deseqResultsDir, findModule, outDir, vargs)

    #
    # create upset plots and find interesection elements
    #
    _step2(logger, outDir, selectedGeneSetsDict, sgc)

    #
    # step 3
    # create data into CIBERSORTx signatureGenes matrix
    #
    logger.info('step 3: create CIBERSORTx signatureGenes matrix')

    csmf = CibersortSignatureMatrixFactory(
        geneSignatureProfilesDataRootDir = sgc.getLocalCachedDir(), 
        groupByGeneCountFilePath = countDataPath,
        colDataFilePath = colDataPath, 
        estimatedScalingFactorsFilePath = esfPath, 
        localCacheDir= outDir, 
        outdir = "ciberSortInput", 
        testSize = None, 
        verbose = False
    )
        
    signatureGeneDF = csmf.getCiberSortSignatueDF()

    logger.info("signatureGeneDF.shape:{}".format(signatureGeneDF.shape)) 
    # self.logger.info( f'signatureGeneDF\n {signatureGeneDF}' )
    pathToSignatureFile = csmf.save()
    logger.warning(f'pathToSignatureFile : {pathToSignatureFile}')

    #
    # step 4
    # create CIBERSORTx mixture matrix
    #
    logger.info('step 4: create CIBERSORTx mixture matrix')

    csmf = CibersortMixtureFactory(
                str(pathToSignatureFile),
                countDataPath,
                colDataPath,
                esfPath,
                outDir
            )
    
    #aedwip prefix = "testCibersortMixtureFactory"
    prefix = None
    # saveDir = os.path.join(outDir, sgc.getLocalCachedDir(), "ciberSortInput")
    saveDir = os.path.join( sgc.getLocalCachedDir(), "ciberSortInput")
    mixturePath, expectedFractionsPath = csmf.saveMixtureAndExpectedFractions(outDir=saveDir, prefixStr=prefix)
    logger.warning(f'mixturePath : {mixturePath}')
    logger.warning(f'expectedFractionsPath : {expectedFractionsPath}')

    randomizedMixtureDF = csmf.randomizeMixture(seed=42)
    randomizedPath = csmf.saveRandomizedMixture(outDir=saveDir, randomizedDF=randomizedMixtureDF, prefixStr=prefix)
    logger.warning(f'randomizedPath : {randomizedPath}')

    logger.warning("END")
    sys.exit(0)

################################################################################
if __name__ == '__main__':
    main()
