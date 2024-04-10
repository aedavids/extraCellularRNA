#
# randomForestHyperparmeterSearch.py
# Andrew E. Davidson
# aedavids@ucsc.edu
# 02/01/2024
#
# ref: extracelluarRNA/intraExtraRNA_POC/jupyterNotebooks/elife/lungCancer/randomForestIntraCellularLungCancerBiomarkersOnExtracellularSamples.ipynb
#

# RandomForestHyperparmeterSearchCLI CommandLine display the doc string 
'''
Random Forest hyperparameter search
'''

import itertools as it
import logging
import numpy as np
import os
import pandas as pd 
from sklearn.ensemble        import RandomForestClassifier
from sklearn.model_selection import  cross_validate
from sklearn.metrics         import recall_score
from sklearn.metrics         import roc_auc_score
from sklearn.metrics         import make_scorer
from sklearn.model_selection import RepeatedStratifiedKFold
import sys

from intraExtraRNA.elifeUtilities import loadElifeTrainingData
from models.randomForestHyperparmeterSearchCLI import RandomForestHyperparmeterSearchCLI

meaningOfLife = 42

# global variables use by PipelineCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-12-07'
__updated__ = '2023-12-07'

################################################################################
def evaluateModel(logger : logging.Logger, 
                  model, 
                  XNP, 
                  yNP, 
                  scoringMetricsDict, 
                  randomSeed=meaningOfLife, ) -> dict:
    '''
    TODO
    '''
    logger.info("BEGIN")

    # define the evaluation procedure
    # crossValidationGenerator = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=meaningOfLife)
    crossValidationGenerator = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=randomSeed)
    
    # evaluate the model and collect the results
    scoresDict = cross_validate(model, XNP, yNP, scoring=scoringMetricsDict, cv=crossValidationGenerator, n_jobs=-1)
    
    logger.info("END")
    return scoresDict

################################################################################
def createScoringMetricsDict(logger : logging.Logger) -> dict:
    '''
    https://stackoverflow.com/a/64540588/4586180
    The specifity is the True Negative Rate which is the same as the True Positive Rate (Recall)
    but for the negative class
    !!! this assume we have a binary classifier.
    
    '''
    logger.info("BEGIN")
    specificity = make_scorer(recall_score, pos_label=0)
    auc = make_scorer( roc_auc_score )

    scoringMetricsDict = {
        'accuracy' : 'accuracy', # TP
        'sensitivity' : 'recall',
        'specificity' : specificity,
        'auc' : auc,

        #'f1' : 'f1'
    } 

    logger.info("END")
    return scoringMetricsDict

################################################################################
def createMaxDepth(logger : logging.Logger,
                          debug :bool=False) -> tuple[ tuple[str, int] ]:
    '''
    TODO
    '''
    logger.info("BEGIN")
    # None is default
    depths = [i for i in range(1,8)] + [None]
    maxDepth = tuple( it.product(['max_depth'], depths) )

    # aediwp
    # modelDict = dict()
    # depths = [i for i in range(1,8)] + [None]
    # for i in depths:
    #     key = str(i)
    #     modelDict[key] = RandomForestClassifier(max_depth=i)
    
    logger.info(f'ret: {maxDepth}')
    logger.info("END")
    return maxDepth

################################################################################
def createMaxFeatures(logger : logging.Logger,
                      XNP: np.array, 
                      bound : int = 5, 
                      debug : bool=False) -> tuple[ tuple[str, int] ]:
    '''
    TODO
    '''
    logger.info("BEGIN")
    # it.product calculate the cartian product
    # example (('max_features', 1), ('max_features', 2), ('max_features', 3))
    nFeatures = XNP.shape[1]
    medianNFeatures = int( nFeatures**0.5 )
    
    if debug :
        bound = 1

    start = medianNFeatures - bound
    if start < 1:
        start = 1
    end = medianNFeatures + bound
    logger.debug(f'AEDWIP start : {start} medianNFeatures: {medianNFeatures} bound : {bound} end : {end}')
    maxFeatures = tuple( it.product( ['max_features'], np.arange(start, end, 1) ) )

    logger.info(f'ret : {maxFeatures}')
    logger.info("END")
    return maxFeatures

################################################################################
def createMaxSamples(logger : logging.Logger,
                     debug : bool=False) -> tuple[ tuple[str, float] ]:
    '''
    TODO
    '''
    logger.info("BEGIN")

    step = 0.1
    start = 0.1
    end = 1.0 + step

    # None is debug
    values = list(np.arange(start, end, step)) + [None]

    if debug:
        end = 0.2 + step
        values = np.arange(start, end, step)

    maxSamples = tuple( it.product( ['max_samples'], values ) )

    logger.info(f'ret : {maxSamples}')
    logger.info('END')
    return maxSamples

################################################################################
def createNumEstimators(logger : logging.Logger,
                              debug :bool=False) ->  tuple[ tuple[str, float] ]:
    '''
    TODO
    '''
    logger.info("BEGIN")

    modelDict = dict()
    
    treeList = [50, 100, 500, 1000, 1500, 2000]
    if debug :
        treeList = [50, 100]

    numTrees = tuple( it.product( ['n_estimators'], treeList ) )

    logger.info(f'ret : {numTrees}')
    logger.info("END")
    return numTrees

################################################################################
def createSearchParameters(logger : logging.Logger,
                           XNP : np.array, 
                           debug : bool) -> list[dict[str, float]]:
    '''
    TODO
    '''
    logger.info("BEGIN")   

    bound = 5
    maxFeatures = createMaxFeatures(logger, XNP, bound, debug )
    maxSamples = createMaxSamples(logger, debug )
    maxTrees = createNumEstimators(logger, debug )
    maxDepth = createMaxDepth(logger, debug )

    #
    # parameters  
    # example parameters items we we only pass maxFeatures, maxSamples to cartian product
    # (('max_features', 1), ('max_samples', 0.1))
    # (('max_features', 1), ('max_samples', 0.5))
    #
    parameters = list( it.product(maxFeatures, maxSamples, maxTrees, maxDepth) )

    kwags = []
    for p in parameters:
        d = dict(p) 
        kwags.append(d)

    logger.info("END")
    return kwags

################################################################################
def tunningFramework(logger : logging.Logger,
                     models : list,
                    XNP : np.array, 
                    yNP : np.array, 
                    scoringMetricsDict : dict) -> dict:
    '''
    TODO
    '''    
    logger.info("BEGIN")
 
    # resultsDict = dict()
    # for key in scoringMetricsDict.keys():
    #     resultsDict[key] = {}

    # create structure to hold all model evaluation results
    resultsDict = dict()
    resultsDict['parameters'] = []
    for metricName in scoringMetricsDict.keys():
        resultsDict[metricName + "_mean"] = list()
        resultsDict[metricName + "_std"] = list()
        
    # for hyperparameterValue, model in maxSampleModelsDict.items():
    for model in models:
        scoresDict = evaluateModel(logger, model, XNP, yNP, scoringMetricsDict)
        resultsDict['parameters'].append( model.aedwipKwags )
        #print(f'\n\n*********scoresDict\n{pp.pformat(scoresDict, indent=4)}\n')
        for metricName in scoringMetricsDict.keys():
            scoreKey = "test_" + metricName
            #print(f'hyperparameterValue: {hyperparameterValue} metricName : {metricName} scoreKey: {scoreKey}')
            values = scoresDict[scoreKey]
            #print(np.round(values, decimals=3))
            mean = np.mean(values)
            std =  np.std(values)
            resultsDict[metricName + "_mean"].append(mean)
            resultsDict[metricName + "_std"].append(std)

            #logger.info(f'{parameters}={hyperparameterValue} {metricName}\tmean : %.3f std : %.3f)' % ( np.mean(values), np.std(values)))
        
        #print(f'scores.keys\n {scores.keys()}')

    logger.info("END")
    return resultsDict

################################################################################
def cleanWhiteSpace(logger,  args = list[str] ):
    ret = []
    for i in range(len(args)):
        c = args[i]
        cc = c.strip()
        logger.info(f'AEDWIP c : xxx{c}xxx cc : xxx{cc}xxx')
        ret.append( cc)

    return ret

################################################################################
def main(inCommandLineArgsList=None):
    '''
    TODO
    '''
    # we only configure logging in main module
    #loglevel = "WARN"
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

    #
    # parse CLI
    #
    cli = RandomForestHyperparmeterSearchCLI(version=__version__ , 
                            author=__author__ ,
                            date=__date__, 
                            update=__updated__)
    
    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    features = cleanWhiteSpace(logger,  cli.args.features )
    outDir   = cli.args.outDir   
    os.makedirs(outDir, exist_ok=True)
    pipelineStageName = cli.args.pipelineStageName
    selectElifeCategories = cleanWhiteSpace(logger, cli.args.elife )

    logger.info(f'features: {features}')
    logger.info(f'outDir: {outDir}')
    logger.info(f'pipelineStageName: {pipelineStageName}')
    logger.info(f'selectElifeCategories: {selectElifeCategories}')

 
    # load training data

    # HUGO_lungGenes, elifeLungGenes, countDF, metaDF, XNP, yNP = loadElifeLungTrainingData(pipelineStageName, features)
    HUGO_lungGenes, elifeLungGenes, countDF, metaDF, XNP, yNP = loadElifeTrainingData(pipelineStageName, 
                                                                                      features, 
                                                                                      selectElifeCategories)

    parameterKwags= createSearchParameters(logger, XNP, debug=False)

    #
    # for each set of search parameters
    # create a model
    #
    models = list()
    for kwags in parameterKwags:
        rm = RandomForestClassifier(**kwags)
        models.append(kwags)

    logger.info(f'len(models): {len(models)}')

    scoringMetricsDict = createScoringMetricsDict(logger) 

    # create structure to hold all model evaluation results
    # resultsDict = dict()
    # resultsDict['parameters'] = list()
    # for metricName in scoringMetricsDict.keys():
    #     resultsDict[metricName + "_mean"] = list()
    #     resultsDict[metricName + "_std"] = list()

    # logger.error(f'********* AEDWIP DEBUG = True')
    searchGridParameters = createSearchParameters(logger, XNP, debug=False)

    #
    # create model for each set of parameters we which to evaluate
    #
    models = list()
    # aedwip = 0
    for kwags in searchGridParameters:
        rf = RandomForestClassifier(**kwags) 
        logger.debug(f'AEDWIP creating rf with  parameters {kwags}')
 

        # make sure we know the models configuration
        rf.aedwipKwags = kwags
        models.append( rf )
               
        # aedwip = aedwip + 1
        # if aedwip == 20:
        #     break

    logger.warning(f'len(models) : {len(models)}')

    resultsDict = tunningFramework(logger, 
                        models,
                       XNP,
                       yNP,
                       scoringMetricsDict)

    df = pd.DataFrame( resultsDict ) 
    # df = df.sort_values(by='sensitivity_mean', ascending=False)
    df = df.sort_values(by='auc_mean', ascending=False)
    # parameters is a dictionary. split into columns
    expandedDictCols = df['parameters'].apply(pd.Series)
    logger.debug(expandedDictCols.to_string())
    # works but puts expand cols at end 
    df = df.join(expandedDictCols)
    df = df.drop('parameters', axis=1)
    #print(df.to_string())


    logger.warning('*************\ntop 10 results sorted by auc_mean')
    logger.warning(df.head(n=10) )

    # logger.error(f'AEDWIP df.info()\n {df.info()}')

    # 
    # If Pandas is writing integer columns as floats when you use the to_csv method,
    # it's likely due to the presence of missing values (NaN) in those columns.
    # In Pandas, NaN is considered a floating-point number, and since a column 
    # can only have one data type, the entire column gets 
    # cast to float if it contains any NaN values.
    # Convert to Integers with Nullable Integer Type
    #
    # this cause problems with read_csv() the col will be type float not in
    # we could put in a marker value like -99999. I do not know what scikitlearn
    # will do?
    #
    # why do we care?
    # some random forest parameters behave very differently is the
    # value is float or an integer
    #
    df['max_depth'] = df['max_depth'].astype('Int32')
    df['max_features'] = df['max_features'].astype('Int32')
    df['n_estimators'] = df['n_estimators'].astype('Int32')

    logger.debug(f'AEDWIP df.info()\n {df.info()}')

    outFile =  f'{outDir}/randomForestHyperparmeterSearch.csv'
    df.to_csv(outFile, index=False)
    logger.warning(f'saved file to {outFile}')

    logger.warning("END")
    sys.exit(0)

################################################################################
def mainCLITest(inCommandLineArgsList=None):
    '''
    hack used to test parsing of vargs
    '''

    # we only configure logging in main module
    # loglevel = "WARN"
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

    cli = RandomForestHyperparmeterSearchCLI(version=__version__ , 
                            author=__author__ ,
                            date=__date__, 
                            update=__updated__)
    
    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    features = cli.args.features 
    outDir   = cli.args.outDir    

    logger.info(f'features: {features}')
    logger.info(f'outDir: {outDir}')

    
################################################################################
if __name__ == '__main__':
    main()
    
    #hack used to test parsing of vargs
    # mainCLITest(
    #     #inCommandLineArgsList=['-h']
    #     inCommandLineArgsList=[
    #         '--features', 'LUAD',
    #         '--outDir', 'tmp'
    #     ]
    # )
