#
# elifeUtilities.py
# Andrew E. Davidson
# aedavids@ucsc.edu
#
'''
public functions 
    loadCounts()
    loadElifeTrainingData()
    loadElifeLungTrainingData()
    loadMetaData()
    searchForMissingMapGenes()
    selectFeatures()
'''
import logging
import pandas as pd
import pprint as pp
import numpy as np
from sklearn.preprocessing import LabelEncoder
 
from analysis.utilities import  findSignatureGenesForPipelineStage
from intraExtraRNA.mapIds import mapHUGO_2_ENSG
from intraExtraRNA.utilities import  selectSamples

logger = logging.getLogger(__name__)

validElifeCategories = {"Colorectal Cancer", "Esophagus Cancer", "Healthy donor", "Liver Cancer", "Lung Cancer", "Stomach Cancer"}

################################################################################
def loadCounts(
        dataRoot : str = "/private/groups/kimlab/alex/data/elife",
        countFile : str = "elife_all_norm_counts_2023-05-18.csv"
        ) -> pd.DataFrame :
    '''
    example of normalized counts
        $ head elife_all_norm_counts_2023-05-18.csv | cut -d , -f 1,2,3,4
        gene,SRR14506659,SRR14506660,SRR14506661
        (A)n,201.6720531316013,110.45077340739788,3722.776395118632
        (AAA)n,0,0,0
        (AAAAAAC)n,0,0,0
        (AAAAAAG)n,0,0,0
        (AAAAAAT)n,0,0,0
        (AAAAAC)n,0,0,0
        (AAAAACA)n,0,0,0
        (AAAAACC)n,0,0,0
        (AAAAACT)n,0,0,0

    returns a data frame with index = sample ids, and columns = genes

    $ pwd
    /private/groups/kimlab/alex/data/elife
    $ cut -d , -f 12 elife_scaled_metaData_2023-05-18.csv| sort | uniq -c
        53 "Colorectal Cancer"
        1 "diagnosis"
        31 "Esophagus Cancer"
        43 "Healthy donor"
        26 "Liver Cancer"
        35 "Lung Cancer"
        36 "Stomach Cancer"    
    '''

    countPath = f'{dataRoot}/{countFile}'
    logger.info(f'countPath : {countPath}')
    countDF = pd.read_csv(countPath, index_col='gene').transpose()
    return countDF


################################################################################
def loadElifeLungTrainingData(
        pipelineStageName : str = "best10CuratedDegree1_ce467ff",
        features : str = "all",
    ) -> tuple[list[str], list[str], pd.DataFrame, pd.DataFrame, np.array, np.array]:
    '''
    Easy of use class. Makes it easy to select lung releated samples from eLife, GTex, and TCGA
    wrapper around loadElifeTrainingData()

    broke this out to maintain random forest hyper parameter search backwards compatiblity
    
    raises: ValueError

    arguments: 
        features: 
            see RandomForestHyperparmeterSearchCLI
            choices=['LUAD', 'LUSC', 'Lung', "all"]

    returns
        (HUGO_lungGenes, elifeLungGenes, countDF, metaDF, XNP, yNP)
    '''
    logger.info("BEGIN")
    if features == "all":# TODO AEDWIP do not hard code
        categories = ['LUAD', 'LUSC', 'Lung']
    else :
        categories = [features]

    logger.debug(f'AEDWIP categories : {categories}')

    selectElifeCategories =  ["Healthy donor", "Lung Cancer"]

    logger.info("END")
    return loadElifeTrainingData(pipelineStageName,
                                 categories,
                                 selectElifeCategories)

################################################################################
def loadElifeTrainingData(
        pipelineStageName : str,
        features : list[str],
        selectElifeCategories : list[str],
    ) -> tuple[list[str], list[str], pd.DataFrame, pd.DataFrame, np.array, np.array]:
    '''

    raises: ValueError

    arguments: 
        features: 
            a list of GTEx or TCGA classes to select features from 

        selectElifeCategories:
            valid items:
                 "Colorectal Cancer", "Esophagus Cancer", "Healthy donor", 
                 "Liver Cancer", "Lung Cancer", "Stomach Cancer"


    returns
        (HUGOGenes, elifeLungGenes, countDF, metaDF, XNP, yNP)
    '''
    logger.info("BEGIN")
    logger.debug(f"AEDWIP pipelineStageName : {pipelineStageName}")
    #logger.error(f'AEDWIP features: {features}')
    countDF = loadCounts()
    metaDF = loadMetaData()

    # if features == "all":# TODO AEDWIP do not hard code
    #     categories = ['LUAD', 'LUSC', 'Lung']
    # else :
    #     categories = [features]

    # logger.debug(f'AEDWIP categories : {categories}')

    # validate elife categories arguments
    errorMsg = None
    for elifeType in selectElifeCategories:
        if not elifeType in validElifeCategories:
            if errorMsg == None:
                errorMsg = "ERROR invalide elife diagnosis: "
            errorMsg = errorMsg + f' xxx{elifeType}xxx expected : {validElifeCategories}'

    if not errorMsg == None:
        logger.error(errorMsg)
        raise ValueError(errorMsg)

    # category="LUAD" 
    # LUADGenes = findSignatureGenesForPipelineStage(category, pipelineStageName, )
    # logger.info(f'LUAD genes:\n{LUADGenes}')

    # category="LUSC"
    # LUSCGenes = findSignatureGenesForPipelineStage(category, pipelineStageName)
    # logger.info(f'\nLUSC genes:\n{LUSCGenes}')

    # # Lung is our healthy control
    # category="Lung"
    # controlGenes = findSignatureGenesForPipelineStage(category, pipelineStageName, )
    # logger.info(f'\n healthy control Lung genes:\n{controlGenes}')

    # HUGO_lungGenes = LUADGenes + LUSCGenes + controlGenes
        
    HUGOGenes = []
    for c in features :
        logger.debug(f'AEDWIP category : {c}')
        genes = findSignatureGenesForPipelineStage(c, pipelineStageName, )
        logger.debug(f'category : {c} genes:\n{genes}')
        HUGOGenes = HUGOGenes + genes

    logger.info(f'len(HUGOGenes) {len(HUGOGenes)}')

    # check for missing biomarkers
    elifeLungGenes, missingElifeGenes = selectFeatures(countDF, HUGOGenes)
    logger.warning( f'len(elifeLungGenes) : {len(elifeLungGenes)}' )
    logger.warning( f'missingElifeGenes\n : {missingElifeGenes}' )

    # for now just drop missing genes
    features = list( set(elifeLungGenes) - set(missingElifeGenes) )
    # features is something like Colon_Sigmoid assert len(features) > 1, f'ERROR: loadElifeTrainingData() {HUGOGenes} do not map to ENSGO '
    # only all has 29 assert len(features) == 29, "ERROR removing missing elife genes"

    # select training data
    #aedwip selectElifeCategories = ["Healthy donor", "Lung Cancer"] aedwip do not hard code , use default value to provide backwards compatiblity
    tmpMetaDF = metaDF.rename( columns={ "diagnosis" : "category"} )
    # display( tmpMetaDF.head() )
    XDF = selectSamples(tmpMetaDF, countDF, selectElifeCategories)

    XDF = XDF.loc[:, features]
    logger.info(f'XDF.shape : {XDF.shape}')

    # create labels
    selectRows = tmpMetaDF.loc[:, 'category'].isin(selectElifeCategories)
    conditionList = tmpMetaDF.loc[selectRows, 'category'].tolist()

    labelEncoder = LabelEncoder()
    yNP = labelEncoder.fit_transform(conditionList)
    XNP = XDF.values

    logger.info("END")
    return (HUGOGenes, elifeLungGenes, countDF, metaDF, XNP, yNP)

################################################################################
def loadMetaData(
                dataRoot : str = "/private/groups/kimlab/alex/data/elife",
                metaFile : str = 'elife_scaled_metaData_2023-05-18.csv',
                columns : list[str] =  ["sample_id", "diagnosis"],
    ) -> pd.DataFrame :
    '''
    example input

    $ head -n 1 elife_scaled_metaData_2023-05-18.csv 
            "sample_id","salmon_frag_length_mean","salmon_frag_length_sd",
            "salmon_num_processed","salmon_num_mapped","salmon_num_decoy_fragments",
            "salmon_num_dovetail_fragments","salmon_num_fragments_filtered_vm",
            "salmon_num_alignments_below_threshold_for_mapped_fragments_vm",
            "salmon_percent_mapped","age","diagnosis","gender","input_vol",
            "dataset","condition"
    '''
    metaDataPath = f'{dataRoot}/{metaFile}'
    metaDF = pd.read_csv(metaDataPath).loc[:, columns]    

    return metaDF

################################################################################
def searchForMissingMapGenes(countDF, genes, refSeq2ENSGDF):
    '''
    elife complete seq count data uses ensembl ids

    We can map the GTEx_TCGA HUGO gene ids to ensembl ids

    The GTEx_TCGA complete seq counts where done using ref gencode.v35

    Some of the gencode.v35 ids do not map to gencode.v39

    linux hack find missing gencode.v35 ENSG00000253339.2 id
    s=/private/groups/kimlab/alex/data/elife/elife_all_norm_counts_2023-05-18.csv
    s.shape : (76556, 225)
    $ grep ENSG00000253339 $s | cut -d , -f 1
    ENSG00000253339.3

    '''
    logger.info("BEGIN")
    cols = countDF.columns

    # make sure all genes exist in elife data
    i = 0
    missingGenesElife = []
    geneSet = set(genes)
    retGeneSet = geneSet.copy()
    for k in  genes:
        if not k in cols:
            i += 1
            logger.info(f'gencode.v35.ucsc.rmsk.tx.to.gene.csv {k} not found in elife gencode.v39.annotation.expanded.tx.to.gene.tsv')
            missingGenesElife.append(k)
    
    # logger.info(f'missing i : {i}')
    if len(missingGenesElife) > 0 :
        selectRows = refSeq2ENSGDF.loc[:, 'ENSG'].isin(missingGenesElife)
        logger.info( f'refSeq2ENSGDF.loc[selectRows, :] :\n{pp.pformat(refSeq2ENSGDF.loc[selectRows, :])}')

        # we can ignore the decimal point. it encode the version number
        # TODO AEDWIP calculate the hackDict
        # use pandas split(expand=true) https://pandas.pydata.org/docs/reference/api/pandas.Series.str.split.html
        hackDict = { # key = v25 value = v39
                     'ENSG00000253339.2' : 'ENSG00000253339.3', # source elife Lung Cancer
                     'ENSG00000267107.8' : 'ENSG00000267107.9', # source elife Lung Cancer

                     'ENSG00000076770.14' : 'ENSG00000076770.15', # Source Colon_Sigmoid elife Colorectal Cancer
                     'ENSG00000111554.14' : 'ENSG00000111554.15', # Source Colon_Sigmoid elife Colorectal Cancer
                     'ENSG00000214944.9' : 'ENSG00000214944.10', # Source Colon_Sigmoid elife Colorectal Cancer

                     'ENSG00000243701.7' : 'ENSG00000243701.8', # COAD elife Colorectal Cancer

                     'ENSG00000137959.16' : 'ENSG00000137959.17', # READ, elife Colorectal Cancer
                     'ENSG00000134202.11' : 'ENSG00000134202.12', # READ, elife Colorectal Cancer

                     #'ENSG00000274031.1' : ' not found ????', # LUSC

                     'ENSG00000225889.9' : 'ENSG00000225889.10', # Stomach
                     'ENSG00000171840.12' : 'ENSG00000171840.13', #Stomach
                     '' : '', # 
                     }
        
        for missing in missingGenesElife:
            if missing in hackDict :
                map2Id = hackDict[missing]
                logger.info(f'mapping {missing} to {map2Id}')
                retGeneSet.add(map2Id)
                retGeneSet.remove(missing)
                missingGenesElife.remove(missing)
    

    logger.info("END")
    return (list(retGeneSet), missingGenesElife)

################################################################################
def selectFeatures(
                    countDF : pd.DataFrame,
                    genes : list[str],
                    txt2GeneFilePath : str =  "/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv"

                   ):
    '''
    TODO
    convert list of HUGO gene ids to ENSG is

    example:
        genes = ['AC011944.1', 'ATP13A4-AS1', 'AC090004.1', 'GXYLT1P3', 'AC126323.6']
    '''
    logger.info("BEGIN")

    geneMapDF = mapHUGO_2_ENSG(txt2GeneFilePath)
    logger.info(f'geneMapDF.shape : {geneMapDF.shape }')

    selectGeneRows = geneMapDF.loc[:, "HUGO"].isin(genes)
    refSeq2ENSGDF = geneMapDF.loc[selectGeneRows.values, ['HUGO', 'ENSG', 'bioType']].drop_duplicates()
    logger.info(f'refSeq2ENSGDF.shape : {refSeq2ENSGDF.shape}')
    logger.info(f'refSeq2ENSGDF.head()\n{refSeq2ENSGDF.head()}')

    elifeGenes = refSeq2ENSGDF.loc[:, "ENSG"].drop_duplicates().tolist()
    logger.info(f'len(elifeGenes) : {len(elifeGenes)}')
    elifeGenes[0:3]

    # check for missing gene loci ids
    missingIdsSet = set(genes) - set(refSeq2ENSGDF.loc[:, 'HUGO'])
    logger.info(f'missingIdsSet\n{missingIdsSet}')

    # make sure missing ids exist in elife count data and add to elifeLungGenes
    elifeGeneSet = set(countDF.columns)
    # missing = set(lungGenes) - set( refSeq2ENSGDF['refSeq'] )
    # logger.info(f'missing :\n{missing}')
    for k in  missingIdsSet:
        if not k in elifeGeneSet:
            logger.info(f'{k} not found')
        else :
            logger.info(f'adding {k} to elifeLungGenes ')    
            elifeGenes.append(k)        

    # assert len(elifeLungGenes) == 30, f'ERROR missing ids expected 30 got {len(elifeLungGenes) }'

    #  Some of the gencode.v35 ids do not map to gencode.v39
    # searchForMissingMapGenes is a hack to correct for mapping issue
    elifeGenes, missingElifeGenes = searchForMissingMapGenes(countDF, elifeGenes, refSeq2ENSGDF)

    logger.info('END')
    return ( elifeGenes, missingElifeGenes )
