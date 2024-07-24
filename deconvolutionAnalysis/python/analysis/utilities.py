#
# utilities.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#

'''
public functions:
    fileNameToDictKey
    findAllCategories
    findAllGenes
    findElementsInIntersectionsWithDegree
    findDir
    findFile
    findIntersectionsWithDegree
    findSetsWithDegree
    findSignatureGenesForPipelineStage
    loadDictionary
    loadlist
    loadPipelineStageIntersectionDict
    saveDict
    saveList
    saveSet
'''

import ast
import fnmatch
import logging
import os
import pandas as pd
import pprint as pp

from pipeline.dataFactory.driver import _countExtraHeaderLines

logger = logging.getLogger(__name__)

################################################################################
def fileNameToDictKey(fileName : str ) ->tuple[str]:
    '''
    returns the degree1 dictionary key for file name

    example filename : 'Whole_Blood_vs_all.results"

    returns ('Whole_Blood', )
    '''
    root = fileName.removesuffix( "_vs_all.results" )
    ret = (root,)

    return ret

################################################################################
def findAllCategories(intersectionDict : dict) -> set :
    '''
    TODO
    '''

    ret = set()
    for intersectionKeys in intersectionDict.keys():
        ret = ret.union( set(intersectionKeys) )

    return ret

################################################################################
def findAllGenes(intersectionDict : dict) -> set :
    '''
    TODO
    '''

    ret = set()
    for intersectionKeys, genes in intersectionDict.items():
        ret = ret.union( set(genes) )

    return ret

################################################################################
def findElementsInIntersectionsWithDegree(
            intersectionDict :dict,
            degree : int) -> set :
    '''
    returns a set containing all the elements from interections with degree
    '''
    
    elements = []
    for intersectionKey, values in intersectionDict.items():
        if len(intersectionKey) == degree:
            elements = elements + list(values)

    return set( elements )

################################################################################
def findDir(rootDir : str, pattern : str) -> list[str]:
    '''
    ref: https://stackoverflow.com/a/1724723/4586180

    arguments:
        pattern: 
            fileName. support for regex
            example: '*.txt'
    '''
    # metricsPath = ! find $path -name "metricsRounded.csv"

    result = []
    for root, dirs, files in os.walk(rootDir):
        for name in dirs:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

################################################################################
def findFile(rootDir : str, pattern : str) -> list[str]:
    '''
    ref: https://stackoverflow.com/a/1724723/4586180

    arguments:
        pattern: 
            fileName. support for regex
            example: '*.txt'
    '''
    # metricsPath = ! find $path -name "metricsRounded.csv"

    result = []
    for root, dirs, files in os.walk(rootDir):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

################################################################################
def findIntersectionsWithDegree(
            intersectionDict :dict,
            degree : int) -> dict :
    '''

    returns a dictionary of all intersections with degree.
    ie a "sub set" of the original dictionary
    '''

    retDict = dict()
    for key, items in intersectionDict.items():
        if len(key) == degree:
            retDict[ key ] = items
    
    return retDict

################################################################################
def findSetsWithDegree(
            intersectionDict :dict,
            degree : int) -> set :
    '''

    returns a set. The elements are the dictionary keys items for keys of length = degree
    '''

    categories = []
    for intersectionKey in intersectionDict.keys():
        if len(intersectionKey) == degree:
            categories = categories + list(intersectionKey)

    return set( categories)

################################################################################
def findSignatureGenesForPipelineStage(
        category : str,
        pipelineStageName : str,         
        #LOCAL_CACHE_DIR : str, 
        deconvolutionRoot : str = "/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category",
        colName : str = "name",
        ascendingDecendingHack : bool = False,
    ) -> list:
    '''
    example 
        category = "LUSC"
        pipelineStageName = "best10CuratedDegree1_ce467ff"
        colName = "name"

    7/2/24 added ascendingDecendingHack. default value = False
    this is a work around for 
    deconvolutionAnalysis/jupyterNotebooks/hyperParameterTunning/ascending-vs.-DescendingBaseMeanGeneSignatureSelection.ipynb
    best10CuratedDegree1 has 2 out dirs best10CuratedDegree1.sh.out-sortOrder-bug and  best10CuratedDegree1.sh.out

    returns a list. The type of the list items is determined by pandas read_csv()
    '''
    # load lung cancer signature, (biomarker), genes
    
    
    #upsetPlotDir = f'{deconvolutionRoot}/{pipelineStageName}/training/best10CuratedDegree1.sh.out/upsetPlot.out'
    searchRoot = f'{deconvolutionRoot}/{pipelineStageName}'
    ciberSortInputList = findDir(searchRoot, "ciberSortInput")

    if not ascendingDecendingHack :
        if len(ciberSortInputList) != 1:
            logger.error(f'pipelineStageName : {pipelineStageName} category : {category} unable to find results in searchRoot: {searchRoot} found {ciberSortInputList}')
            return 

    # find the deseq results file used to create the signature matrix
    deseqResultsDir = os.path.dirname( ciberSortInputList[0] )
    resultFile = os.path.join(deseqResultsDir, category + '_vs_all.results')
    logger.info(f'resultFile : {resultFile}')
    numRowsToSkip = _countExtraHeaderLines(resultFile)
    logger.debug(f'numRowsToSkip : {numRowsToSkip}')

    deseqDF = pd.read_csv(resultFile, skiprows=numRowsToSkip)
    retList = deseqDF.loc[:, colName].tolist()

    return retList

################################################################################
def loadDictionary(intersectionDictPath : str) -> dict:
    '''
    TODO
    '''
    with open(intersectionDictPath) as f: 
            data = f.read() 

    # we can not convert the intersectionData to a DataFrame
    # the size of the interesections is not constant. 
    # we would need to pad the dictionary values
    retDict = ast.literal_eval(data)

    return retDict

################################################################################
def loadList(listPath : str) -> list:
    '''
    TODO does not work with SaveList( isSingleItemLine=True )
    '''
    with open(listPath, "r") as f: 
        data = f.read() 

    # we can not convert the intersectionData to a DataFrame
    # the size of the interesections is not constant. 
    # we would need to pad the dictionary values
    ret = ast.literal_eval(data)

    return data

################################################################################
def loadPipelineStageIntersectionDict(
        pipelineStageName : str,         
        #LOCAL_CACHE_DIR : str, 
        deconvolutionRoot : str = "/private/groups/kimlab/aedavids/deconvolution/1vsAll-~gender_category",
    ) -> dict[tuple[str], list[str]]:
    '''
    example 
        pipelineStageName = "best10CuratedDegree1_ce467ff"
    '''
    # load lung cancer signature, (biomarker), genes
    
    
    #upsetPlotDir = f'{deconvolutionRoot}/{pipelineStageName}/training/best10CuratedDegree1.sh.out/upsetPlot.out'
    searchRoot = f'{deconvolutionRoot}/{pipelineStageName}'
    upsetPlotDirList = findDir(searchRoot, "upsetPlot.out")
    if len(upsetPlotDirList) != 1:
        logger.error(f'unable to find upsetPlot.out in {searchRoot} found {upsetPlotDirList}')
        return 

    upsetPlotDir = upsetPlotDirList[0]
    # the deconvolution gene signature matrix gene name are the values in the intersection dictionary
    #intersectionDictPath = f'{upsetPlotDir}/best10_from_best500FindAllDegree1_wl500.intersection.dict'
    intersectionDictPathList = findFile(upsetPlotDir, "*intersection.dict")
    if len(intersectionDictPathList) != 1:
            logger.error(f'unable to find "*intersection.dict" in {upsetPlotDir} found {intersectionDictPathList}')
            return 

    intersectionDictPath = intersectionDictPathList[0]
    # from intraExtraRNA.utilities import load,
    #intersectionDictPath = load(intersectionDictPath, localCacheDir=LOCAL_CACHE_DIR, verbose=True)

    #print(intersectionDictPath)
    intersectionDict = loadDictionary(intersectionDictPath)

    return intersectionDict

################################################################################
def loadSet(setPath : str) -> set:
    '''
    TODO
    '''
    with open(setPath, "r") as f: 
        data = f.read() 

    # we can not convert the intersectionData to a DataFrame
    # the size of the interesections is not constant. 
    # we would need to pad the dictionary values
    ret = ast.literal_eval(data)

    return ret

################################################################################
def saveDictionary(dictPath : str, d: dict):
    '''
    TODO
    '''
    with open(dictPath, "w") as f: 
        f.write(pp.pformat(d,  indent=4, sort_dicts=True)) 
        f.write("\n")

################################################################################
def saveList(listPath : str, l : list, isSingleItemLine : bool=False):
    '''
    TODO
    '''
    with open(listPath, "w") as f:
        if not isSingleItemLine:
            f.write(str(l))
        else:
            for value in l:
                f.write(str(value) +"\n")
    
################################################################################
def saveSet(setPath : str, s: set):
    '''
    TODO
    '''
    with open(setPath, "w") as f: 
        f.write(str(s)) 

