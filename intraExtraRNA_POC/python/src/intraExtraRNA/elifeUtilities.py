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
    t = loadElifeTrainingData(pipelineStageName,
                                 categories,
                                 selectElifeCategories)

    HUGOGenes, elifeLungGenes, missingGenes, countDF, metaDF, XNP, yNP, labelEncoder, mapDF = t

    # for backwards compatibility do not return missing Genes
    return (HUGOGenes, elifeLungGenes, countDF, metaDF, XNP, yNP, labelEncoder, mapDF)

################################################################################
def loadElifeTrainingData(
        pipelineStageName : str,
        features : list[str],
        selectElifeCategories : list[str],
    ) -> tuple[list[str], list[str], list[str], pd.DataFrame, pd.DataFrame, np.array, np.array, LabelEncoder, pd.DataFrame]:
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
        (HUGOGenes, elifeLungGenes, missingGenes, countDF, metaDF, XNP, yNP, labelEncoder, mapDF)

        mapDF shows how v35 HUGO names where mapped to v39 ENSG names
    '''
    logger.info("BEGIN")
    logger.debug(f"AEDWIP pipelineStageName : {pipelineStageName}")
    #logger.error(f'AEDWIP features: {features}')
    transposedCountsDF = loadCounts()
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
    elifeLungGenes, missingElifeGenes, mapDF = selectFeatures(transposedCountsDF, HUGOGenes)
    logger.info( f'len(elifeLungGenes) : {len(elifeLungGenes)}' )
    if len(missingElifeGenes) > 0:
        logger.warning( f'missingElifeGenes\n : {missingElifeGenes}' )

    # for now just drop missing genes
    features = list( set(elifeLungGenes) - set(missingElifeGenes) )
    # features is something like Colon_Sigmoid assert len(features) > 1, f'ERROR: loadElifeTrainingData() {HUGOGenes} do not map to ENSGO '
    # only all has 29 assert len(features) == 29, "ERROR removing missing elife genes"

    # select training data
    #aedwip selectElifeCategories = ["Healthy donor", "Lung Cancer"] aedwip do not hard code , use default value to provide backwards compatiblity
    tmpMetaDF = metaDF.rename( columns={ "diagnosis" : "category"} )
    # display( tmpMetaDF.head() )
    XDF = selectSamples(tmpMetaDF, transposedCountsDF, selectElifeCategories)

    XDF = XDF.loc[:, features]
    logger.info(f'XDF.shape : {XDF.shape}')

    # create labels
    selectRows = tmpMetaDF.loc[:, 'category'].isin(selectElifeCategories)
    conditionList = tmpMetaDF.loc[selectRows, 'category'].tolist()

    labelEncoder = LabelEncoder()
    yNP = labelEncoder.fit_transform(conditionList)
    XNP = XDF.values

    logger.info("END")
    return (HUGOGenes, elifeLungGenes, missingElifeGenes, transposedCountsDF, metaDF, XNP, yNP, labelEncoder, mapDF)

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
# def DEPRECATED_v35biomarkerGeneMapDF(countDF, genes, refSeq2ENSGDF):
#     '''
#     elife complete seq count data uses ensembl ids

#     We can map the GTEx_TCGA HUGO gene ids to ensembl ids

#     The GTEx_TCGA complete seq counts where done using ref gencode.v35

#     Some of the gencode.v35 ids do not map to gencode.v39

#     linux hack find missing gencode.v35 ENSG00000253339.2 id
#     s=/private/groups/kimlab/alex/data/elife/elife_all_norm_counts_2023-05-18.csv
#     s.shape : (76556, 225)
#     $ grep ENSG00000253339 $s | cut -d , -f 1
#     ENSG00000253339.3

#     '''
#     logger.info("BEGIN")
#     cols = countDF.columns

#     # make sure all genes exist in elife data
#     i = 0
#     missingGenesElife = []
#     geneSet = set(genes)
#     retGeneSet = geneSet.copy()
#     for k in  genes:
#         if not k in cols:
#             i += 1
#             logger.info(f'gencode.v35.ucsc.rmsk.tx.to.gene.csv {k} not found in elife gencode.v39.annotation.expanded.tx.to.gene.tsv')
#             missingGenesElife.append(k)
    
#     # logger.info(f'missing i : {i}')
#     if len(missingGenesElife) > 0 :
#         selectRows = refSeq2ENSGDF.loc[:, 'ENSG'].isin(missingGenesElife)
#         logger.info( f'refSeq2ENSGDF.loc[selectRows, :] :\n{pp.pformat(refSeq2ENSGDF.loc[selectRows, :])}')

#         # we can ignore the decimal point. it encode the version number
#         # TODO AEDWIP calculate the hackDict
#         # use pandas split(expand=true) https://pandas.pydata.org/docs/reference/api/pandas.Series.str.split.html
#         hackDict = { # key = v25 value = v39
#                      'ENSG00000253339.2' : 'ENSG00000253339.3', # source elife Lung Cancer
#                      'ENSG00000267107.8' : 'ENSG00000267107.9', # source elife Lung Cancer

#                      'ENSG00000076770.14' : 'ENSG00000076770.15', # Source Colon_Sigmoid elife Colorectal Cancer
#                      'ENSG00000111554.14' : 'ENSG00000111554.15', # Source Colon_Sigmoid elife Colorectal Cancer
#                      'ENSG00000214944.9' : 'ENSG00000214944.10', # Source Colon_Sigmoid elife Colorectal Cancer

#                      'ENSG00000243701.7' : 'ENSG00000243701.8', # COAD elife Colorectal Cancer

#                      'ENSG00000137959.16' : 'ENSG00000137959.17', # READ, elife Colorectal Cancer
#                      'ENSG00000134202.11' : 'ENSG00000134202.12', # READ, elife Colorectal Cancer

#                      #'ENSG00000274031.1' : ' not found ????', # LUSC

#                      'ENSG00000225889.9' : 'ENSG00000225889.10', # Stomach
#                      'ENSG00000171840.12' : 'ENSG00000171840.13', #Stomach
#                      '' : '', # 
#                      }
        
#         for missing in missingGenesElife:
#             if missing in hackDict :
#                 map2Id = hackDict[missing]
#                 retGeneSet.add(map2Id)
#                 retGeneSet.remove(missing)
#                 missingGenesElife.remove(missing)
    

#     logger.info("END")
#     return (list(retGeneSet), missingGenesElife)

################################################################################
def searchForMissingMapGenes(countDF : pd.DataFrame, 
                             v35Biomarker_Genes : list[str], 
                             v35RefBiomarkersDF : pd.DataFrame,):
    '''
    elife complete seq count data uses ENSGO ids

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
    # countDF has the gene ids we need to map to
    cols = countDF.columns

    # make sure all genes exist in elife data
    i = 0
    missingGenesElife = []
    v35Biomarker_GenesSet = set(v35Biomarker_Genes)
    retGeneSet = v35Biomarker_GenesSet.copy()

    # find v35 genes that do not map directly into v39
    # typically it is because the ENSGO version number has changed
    for k in  v35Biomarker_GenesSet:
        if not k in cols:
            i += 1
            logger.info(f'gencode.v35.ucsc.rmsk.tx.to.gene.csv {k} not found in elife gencode.v39.annotation.expanded.tx.to.gene.tsv')
            missingGenesElife.append(k)
    
    # logger.info(f'missing i : {i}')
    if len(missingGenesElife) > 0 :
        logger.info(f'v35 genes that did not map directly to v39: {missingGenesElife}')
        selectRows = v35RefBiomarkersDF.loc[:, 'ENSG'].isin(missingGenesElife)
        logger.info( f'v35RefBiomarkersDF.loc[selectRows, :] :\n{pp.pformat(v35RefBiomarkersDF.loc[selectRows, :])}')

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

                    # hack for 
                    # /private/groups/kimlab/vikas/nanopore/promethion/barretts/analysis/complete-seq/R_results//normalized_counts_for_andy.csv
                    # see /private/groups/kimlab/aedavids/londonCalling2024/data/mapHack.sh
                    'ENSG00000169299.14' : 'ENSG00000169299.14',
                    'ENSG00000051596.10' : 'ENSG00000051596.10',
                    'ENSG00000157379.14' : 'ENSG00000157379.14',
                    'ENSG00000165914.15' : 'ENSG00000165914.15',
                    'ENSG00000129235.11' : 'ENSG00000129235.11',
                    'ENSG00000160703.16' : 'ENSG00000160703.16',
                    'ENSG00000148429.14' : 'ENSG00000148429.15',
                    'ENSG00000129219.14' : 'ENSG00000129219.14',
                    'ENSG00000180667.10' : 'ENSG00000180667.11',
                    'ENSG00000082269.16' : 'ENSG00000082269.17',

                    # hack for ESCA
                    'ENSG00000217325.2' : 'ENSG00000217325.2',
                    'ENSG00000261739.2' : 'ENSG00000261739.2',
                    'ENSG00000224126.2' : 'ENSG00000224126.2',
                    'ENSG00000267125.2' : 'ENSG00000267125.2',
                    'ENSG00000268120.1' : 'ENSG00000268120.1',
                    'ENSG00000203952.9' : 'ENSG00000203952.9',

                    # hack for nanopore ESCA
                    'LTR106' : 'LTR106_Mam'
                     }

        logger.info(f'AEDWIP BEING are there dups? missingGenesElife\n{missingGenesElife}')
        logger.info(f'AEDWIP hackDict:\n{pp.pformat(hackDict, indent=4)}')
        logger.info(f'AEDWIP BEING len(missingGenesElife) : {len(missingGenesElife)}')
        logger.info(f'AEDWIP BEING len(retGeneSet) : {len(retGeneSet)}')
        logger.info(f'AEDWIP BEING before hackDict correction retGeneSet  : {pp.pformat(retGeneSet, indent=4) }')

        i = 0
        tmpGenes = missingGenesElife.copy()
        for missing in tmpGenes:
            missing = missing.strip()
            logger.info(f'AEDWIP i: {i} test missing {missing} is in hack')
            i = i + 1
            if missing in hackDict :
                map2Id = hackDict[missing]
                map2Id = map2Id.strip()
                logger.info(f'mapping {missing} to {map2Id}')
                if missing != map2Id:
                    logger.info(f'{missing} != {map2Id}')
                    retGeneSet.add(map2Id)
                    retGeneSet.remove(missing)
                missingGenesElife.remove(missing)
            else:
                logger.info(f'AEDWIP unable find mapping for {missing}')
                retGeneSet.remove(missing)
    
        logger.info(f'AEDWIP END missingGenesElife\n{missingGenesElife}')
        logger.info(f'AEDWIP END retGeneSet\n{retGeneSet}')
        logger.info(f'AEDWIP END list(retGeneSet)\n{list(retGeneSet)}')


    logger.info("END")
    return (list(retGeneSet), missingGenesElife)

################################################################################
def selectFeatures(
                    countDF : pd.DataFrame,
                    genes : list[str],
                    txt2GeneFilePath : str =  "/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv"
                   ) -> tuple[list[str], list[str], pd.DataFrame]:
    '''
    TODO
    convert list of HUGO gene ids to ENSG is

    example:
        genes = ['AC011944.1', 'ATP13A4-AS1', 'AC090004.1', 'GXYLT1P3', 'AC126323.6']

    returns 
     ( elifeGenes, missingElifeGenes, retMapDF )
    '''
    logger.info("BEGIN")
    logger.info(f'genes: {genes}')

    #
    # load v35 mapping file and select genes of interest
    #
    # mapHUGO_2_ENSG returns a data frame with columns [HUGO ENSG bioType base version]
    # the data frame returns has denormalized entries for repeat gene name
    #  repeats do not have ENSG name we just use what ever name was provide in the 
    # rmask portion of the mapping file
    #
    # denormalization allows use to search for HUGO, ENSG and repeats using the same data frame
    # see mapHUGO_2_ENSG doc string for details
    v35GeneMapDF = mapHUGO_2_ENSG(txt2GeneFilePath)
    logger.info(f'v35GeneMapDF.shape : {v35GeneMapDF.shape }')
    selectGeneRows = v35GeneMapDF.loc[:, "HUGO"].isin(genes)
    logger.info(f'sum(selectGeneRows): {sum(selectGeneRows)}')
    v35RefBiomarkersDF = v35GeneMapDF.loc[selectGeneRows, :]
    logger.info(f'v35RefBiomarkersDF.shape : {v35RefBiomarkersDF.shape}')

    # elife transcripts counts where create using v39 repeat mask
    # load v39 mapping file and select genes of interest
    # i.e the v35 ENSG base == v39 ENSG base
    logger.error(f'AEDWIP do not hard code ')
    aedwip = "/private/groups/kimlab/genomes.annotations/genomes.annotations/gencode.39/gencode.v39.ucsc.rmsk.tx.to.gene.csv"
    v39GeneMapDF= mapHUGO_2_ENSG(aedwip)
    logger.info(f'v39GeneMapDF.shape : {v39GeneMapDF.shape}')

    v35ENSGBaseGeneName = v35RefBiomarkersDF.loc[:,"base"]
    logger.info(f'v35ENSGBaseGeneName + repeats.shape: {v35ENSGBaseGeneName.shape}')
    selectV39Rows = v39GeneMapDF.loc[:,"base"].isin( v35ENSGBaseGeneName)
    v39RefBiomarkersDF = v39GeneMapDF.loc[selectV39Rows, :]
    logger.info(f'v39RefBiomarkersDF.shape : {v39RefBiomarkersDF.shape}')


    # Because of biotype, starting with a list of 10 genes
    # v35RefBiomarkersDF and v39RefBiomarkersDF with have shapes around (40k, 5)
    # join takes for ever. drop duplicate caused by biotypes to reduce size
    v35MapDF = v35RefBiomarkersDF.loc[:, ['HUGO', 'ENSG', 'base', 'version'] ].drop_duplicates()
    logger.info(f'v35MapDF.shape : {v35MapDF.shape}')
    v39MapDF = v39RefBiomarkersDF.loc[:, ['HUGO', 'ENSG', 'base', 'version'] ].drop_duplicates()
    logger.info(f'v39MapDF.shape : {v39MapDF.shape}')

    #mapDF = pd.merge( v35RefBiomarkersDF, v39RefBiomarkersDF, how='inner', on="base", suffixes=('_v35', '_v39') )
    mapDF = pd.merge( v35MapDF, v39MapDF, how='inner', on="base", suffixes=('_v35', '_v39') )

    # get a list of input genes in HUGO mapped to ENSG 
    elifeGenes = list( mapDF.loc[:, 'ENSG_v39'].drop_duplicates() )

    # search for unmapped, i.e missing genes case 1:  gene is not in ENSG_v35
    missingInV35Set = set(genes) - set(v35RefBiomarkersDF.loc[:,'HUGO'])
    if len(missingInV35Set) > 0 :
        logger.warning(f"missingInV35Set: {missingInV35Set}")

    # search for unmapped, i.e missing genes case 2 : gene in ENSG_v35 not in ENSG_v39
    missingInV39Set = set(v35RefBiomarkersDF.loc[:,'base']) - set(v39RefBiomarkersDF.loc[:,'base'])
    if len(missingInV39Set) > 0 :
        logger.warning(f"missingInV39Set: {missingInV39Set}")

    missingElifeGenes = list( missingInV35Set.union(missingInV39Set) )


    # #v35RefBiomarkersDF = v35biomarkerGeneMapDF.loc[selectGeneRows.values, ['HUGO', 'ENSG', 'bioType']].drop_duplicates()
    # # drop duplicate rows
    # #v35RefBiomarkersDF = v35biomarkerGeneMapDF.loc[selectGeneRows.values, ['HUGO', 'ENSG', 'bioType', 'base', 'version']].drop_duplicates()
    # logger.info(f'v35RefBiomarkersDF.shape : {v35RefBiomarkersDF.shape}')
    # #logger.info(f'v35RefBiomarkersDF\n{v35RefBiomarkersDF}')


    # v35ENSG_BiomarkerGeneNames = v35RefBiomarkersDF.loc[:, "ENSG"].drop_duplicates().tolist()
    # logger.info(f'len(v35ENSG_BiomarkerGeneNames) : {len(v35ENSG_BiomarkerGeneNames)}')
    # logger.info(f'v35ENSG_BiomarkerGeneNames : {v35ENSG_BiomarkerGeneNames}')

    # # check for missing gene loci ids
    # missingV35HUGONamesSet = set(genes) - set(v35RefBiomarkersDF.loc[:, 'HUGO'])
    # logger.warning(f'missingV35HUGO HUGO or repeat ids from {txt2GeneFilePath}\n{missingV35HUGONamesSet}')


    # logger.error(f'AEDWIP do not hard code ')
    # aedwip = "/private/groups/kimlab/genomes.annotations/genomes.annotations/gencode.39/gencode.v39.ucsc.rmsk.tx.to.gene.csv"
    # v39RefBiomarkersDF = mapHUGO_2_ENSG(aedwip)
    # logger.info(f'v39RefBiomarkersDF.shape : {v39RefBiomarkersDF.shape}')


    # # make sure missing ids exist in elife count data and add to elifeLungGenes
    # allElifeGeneSet = set(countDF.columns)
    # # missing = set(lungGenes) - set( refSeq2ENSGDF['refSeq'] )
    # # logger.info(f'missing :\n{missing}')
    # for k in  missingV35IdsSet:
    #     if not k in allElifeGeneSet:
    #         logger.info(f'{k} not found')
    #     else :
    #         logger.info(f'adding {k} to elifeLungGenes ')    
    #         v35Biomarker_Genes.append(k)        

    # # assert len(elifeLungGenes) == 30, f'ERROR missing ids expected 30 got {len(elifeLungGenes) }'

    # #  Some of the gencode.v35 ids do not map to gencode.v39
    # # searchForMissingMapGenes is a hack to correct for mapping issue
    # elifeGenes, missingElifeGenes = searchForMissingMapGenes(countDF, v35Biomarker_Genes, v35RefBiomarkersDF)

    retMapDF = mapDF.loc[:, ["HUGO_v35", "ENSG_v35", "ENSG_v39"]]
    logger.info('END')
    return ( elifeGenes, missingElifeGenes, retMapDF )



#############
############# searchForMissingMapGenes prototype
###############

# if we can find a file like /private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv
# for gencode v39 we could expand the ENSG and join on base. most mismapps are due to change in version number
# does gencode.v39.annotation.expanded.tx.to.gene.tsv exist?

# genes
# ['ENSG00000214944.9', 'ENSG00000111554.14', 'ENSG00000076770.14']
# missingGenesDF = pd.DataFrame({"ENSG":genes})

# missingGenesDF
#                  ENSG
# 0   ENSG00000214944.9
# 1  ENSG00000111554.14
# 2  ENSG00000076770.14
# missingGenesDF[["base", "version"]] = missingGenesDF.loc[:,"ENSG"].str.split(".", expand=True)



# missingGenesDF
#                  ENSG             base version
# 0   ENSG00000214944.9  ENSG00000214944       9
# 1  ENSG00000111554.14  ENSG00000111554      14
# 2  ENSG00000076770.14  ENSG00000076770      14



#  refSeq2ENSGDF.head()
#             HUGO                ENSG                  bioType
# 64075   ARHGEF28   ENSG00000214944.9           protein_coding
# 64081   ARHGEF28   ENSG00000214944.9          retained_intron
# 64083   ARHGEF28   ENSG00000214944.9     processed_transcript
# 139028      MDM1  ENSG00000111554.14  nonsense_mediated_decay
# 139030      MDM1  ENSG00000111554.14           protein_coding





# refSeq2ENSGDF[["base", "version"]] = refSeq2ENSGDF.loc[:,"ENSG"].str.split(".", expand=True)

# refSeq2ENSGDF.head()
#             HUGO                ENSG                  bioType             base version
# 64075   ARHGEF28   ENSG00000214944.9           protein_coding  ENSG00000214944       9
# 64081   ARHGEF28   ENSG00000214944.9          retained_intron  ENSG00000214944       9
# 64083   ARHGEF28   ENSG00000214944.9     processed_transcript  ENSG00000214944       9
# 139028      MDM1  ENSG00000111554.14  nonsense_mediated_decay  ENSG00000111554      14
# 139030      MDM1  ENSG00000111554.14           protein_coding  ENSG00000111554      14



# xxxDF = pd.merge

# xxxDF = pd.merge(missingGenesDF,  refSeq2ENSGDF, on='base', how='inner')


# xxx
# Traceback (most recent call last):
#   File "<string>", line 1, in <module>
# NameError: name 'xxx' is not defined
# xxxDF
#                ENSG_x             base version_x      HUGO              ENSG_y                  bioType version_y
# 0   ENSG00000214944.9  ENSG00000214944         9  ARHGEF28   ENSG00000214944.9           protein_coding         9
# 1   ENSG00000214944.9  ENSG00000214944         9  ARHGEF28   ENSG00000214944.9          retained_intron         9
# 2   ENSG00000214944.9  ENSG00000214944         9  ARHGEF28   ENSG00000214944.9     processed_transcript         9
# 3  ENSG00000111554.14  ENSG00000111554        14      MDM1  ENSG00000111554.14  nonsense_mediated_decay        14
# 4  ENSG00000111554.14  ENSG00000111554        14      MDM1  ENSG00000111554.14           protein_coding        14
# 5  ENSG00000111554.14  ENSG00000111554        14      MDM1  ENSG00000111554.14          retained_intron        14
# 6  ENSG00000111554.14  ENSG00000111554        14      MDM1  ENSG00000111554.14     processed_transcript        14
# 7  ENSG00000076770.14  ENSG00000076770        14     MBNL3  ENSG00000076770.14           protein_coding        14
# 8  ENSG00000076770.14  ENSG00000076770        14     MBNL3  ENSG00000076770.14     processed_transcript        14


# missingGenesElife
# ['ENSG00000214944.9']

# missingGenesElife
# ['ENSG00000214944.9', 'ENSG00000111554.14', 'ENSG00000076770.14']

# xxxDF
#                ENSG_x             base version_x      HUGO              ENSG_y                  bioType version_y
# 0   ENSG00000214944.9  ENSG00000214944         9  ARHGEF28   ENSG00000214944.9           protein_coding         9
# 1   ENSG00000214944.9  ENSG00000214944         9  ARHGEF28   ENSG00000214944.9          retained_intron         9
# 2   ENSG00000214944.9  ENSG00000214944         9  ARHGEF28   ENSG00000214944.9     processed_transcript         9
# 3  ENSG00000111554.14  ENSG00000111554        14      MDM1  ENSG00000111554.14  nonsense_mediated_decay        14
# 4  ENSG00000111554.14  ENSG00000111554        14      MDM1  ENSG00000111554.14           protein_coding        14
# 5  ENSG00000111554.14  ENSG00000111554        14      MDM1  ENSG00000111554.14          retained_intron        14
# 6  ENSG00000111554.14  ENSG00000111554        14      MDM1  ENSG00000111554.14     processed_transcript        14
# 7  ENSG00000076770.14  ENSG00000076770        14     MBNL3  ENSG00000076770.14           protein_coding        14
# 8  ENSG00000076770.14  ENSG00000076770        14     MBNL3  ENSG00000076770.14     processed_transcript        14
