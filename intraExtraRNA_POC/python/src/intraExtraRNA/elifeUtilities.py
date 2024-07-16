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
def fixBest10CuratedDegree1_ce467ff(
        elifeLungGenes : list[str], 
        missingElifeGenes : list[str] ) -> tuple[ list[str], list[str] ]:
    '''
    best10CuratedDegree1_ce467f has 716 genes. for unknow reason our current v35 to v39 mapping
    has trouble with a few of these genes. We can not make predictions on elife sample that
    are missing features our random forest was trained on

    ref: intraExtraRNA_POC/jupyterNotebooks/elife/randomForestGTEx_TCGAGeneSignature.ipynb
    
    missing genes we correct for
    ['ENSG00000263264', 'ENSG00000244693', 'ENSG00000288380', 'ENSG00000274031', 'ENSG00000182378.15_PAR_Y']


    returns adjusted version of elifeLungGenes and missingElifeGenes

    '''
    # missingGenes:
    # ['ENSG00000263264', 'ENSG00000244693', 'ENSG00000288380', 'ENSG00000274031', 'ENSG00000182378.15_PAR_Y']

    # (extraCellularRNA) aedavids@mustard $ pwd
    # /private/groups/kimlab/genomes.annotations/gencode.35

    #  $ echo $grepPattern 
    # ENSG00000263264\|ENSG00000244693\|ENSG00000288380\|ENSG00000274031\|ENSG00000182378

    #  grep $grepPattern gencode.v35.ucsc.rmsk.tx.to.gene.csv

    # 'ENSG00000263264'
    # ENST00000576789.1|ENSG00000263264.2|OTTHUMG00000177340.5|OTTHUMT00000436334.2|AC119396.1-201|AC119396.1|1558|protein_coding|,AC119396.1
    # ENST00000671891.1|ENSG00000263264.2|OTTHUMG00000177340.5|OTTHUMT00000529342.1|AC119396.1-202|AC119396.1|1347|protein_coding|,AC119396.1

    # 'ENSG00000244693'
    # ENST00000487179.1|ENSG00000244693.1|OTTHUMG00000158009.1|OTTHUMT00000349996.1|CTAGE8-201|CTAGE8|2588|protein_coding|,CTAGE8

    # 'ENSG00000288380'
    # ENST00000324803.6|ENSG00000288380.1|OTTHUMG00000195956.1|-|AC118281.1-201|AC118281.1|4441|protein_coding|,AC118281.1

    # 'ENSG00000274031'
    # ENST00000618027.1|ENSG00000274031.1|OTTHUMG00000186939.1|OTTHUMT00000474055.1|AC092140.2-201|AC092140.2|521|lncRNA|,AC092140.2


    # 'ENSG00000182378.15_PAR_Y'
    # mapDF
    # idx   HUGO_v35	ENSG_v35	ENSG_v39
    # 670	PLCXD1	ENSG00000182378.14	ENSG00000182378.15
    # 671	PLCXD1	ENSG00000182378.14	ENSG00000182378.15_PAR_Y

    mappingBugDict = {
        'ENSG00000263264' : 'ENSG00000263264.2',
        "ENSG00000244693" : "ENSG00000269693.1", 
        'ENSG00000288380' : 'ENSG00000288380.1',
        'ENSG00000274031' : 'ENSG00000274031.1',
        "ENSG00000182378.15_PAR_Y" : "ENSG00000182378.15",
    }

    for k in mappingBugDict.keys():
        fix = mappingBugDict[k]

        if k in elifeLungGenes:
            logger.warning(f'elifeLungGenes replacing {k} with {fix} elifeLungGenes')
            elifeLungGenes.remove(k) 
            elifeLungGenes.append(fix)         

        if k in missingElifeGenes:
            missingElifeGenes.remove( k )  
            elifeLungGenes.append(fix)  

    # aedwip
    # AEDWIP_BUG_GENE = 'ENSG00000182378.15_PAR_Y'
    # AEDWIP_BUG_GENE_FIX = 'ENSG00000182378.15'
    # if  AEDWIP_BUG_GENE in elifeLungGenes:
    #     # mapDF
    #     # idx   HUGO_v35	ENSG_v35	ENSG_v39
    #     # 670	PLCXD1	ENSG00000182378.14	ENSG00000182378.15
    #     # 671	PLCXD1	ENSG00000182378.14	ENSG00000182378.15_PAR_Y

    #     aedwip 589	AC010422.6	ENSG00000269693.1	ENSG00000269693.1
    #     aedwip 124	LMOD3	ENSG00000163380.16	ENSG00000163380.16

    #     logger.warning(f'elifeLungGenes replacing {AEDWIP_BUG_GENE} with {AEDWIP_BUG_GENE_FIX} elifeLungGenes')
    #     elifeLungGenes.remove(AEDWIP_BUG_GENE) 
    #     elifeLungGenes.append(AEDWIP_BUG_GENE_FIX)
    #     #missingElifeGenes.append( AEDWIP_BUG_GENE )

    return (elifeLungGenes, missingElifeGenes)

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
        (HUGOGenes, elifeGenes, missingGenes, countDF, metaDF, XNP, yNP, labelEncoder, mapDF)

        The columns of XNP are same as the elifeGenes list

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
    elifeLungGenes, missingElifeGenes, mapDF = selectFeatures( HUGOGenes)

    elifeLungGenes, missingElifeGenes = fixBest10CuratedDegree1_ce467ff(elifeLungGenes, missingElifeGenes )

    logger.info( f'len(elifeLungGenes) : {len(elifeLungGenes)}' )
    if len(missingElifeGenes) > 0:
        logger.warning( f'missingElifeGenes\n : {missingElifeGenes}' )

    # for now just drop missing genes
    #features = list( set(elifeLungGenes) - set(missingElifeGenes) )
    # features is something like Colon_Sigmoid assert len(features) > 1, f'ERROR: loadElifeTrainingData() {HUGOGenes} do not map to ENSGO '
    # only all has 29 assert len(features) == 29, "ERROR removing missing elife genes"

    # select training data
    #aedwip selectElifeCategories = ["Healthy donor", "Lung Cancer"] aedwip do not hard code , use default value to provide backwards compatiblity
    tmpMetaDF = metaDF.rename( columns={ "diagnosis" : "category"} )
    # display( tmpMetaDF.head() )
    XDF = selectSamples(tmpMetaDF, transposedCountsDF, selectElifeCategories)
    logger.info(f'AEDWIP XDF.shape : {XDF.shape}')
    logger.info(f'AEDWIP elifeLungGenes : \n{elifeLungGenes}')
    # 77k columns do not log logger.info(f'AEDWIP XDF.columns : \n{XDF.columns}')

 

    # XDF = XDF.loc[:, features]
    XDF = XDF.loc[:, elifeLungGenes]
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
def selectFeatures(
                    genes : list[str],
                    oldTxt2GeneFilePath : str =  "/private/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv",
                    newTxt2GeneFilePath : str =  "/private/groups/kimlab/genomes.annotations/genomes.annotations/gencode.39/gencode.v39.ucsc.rmsk.tx.to.gene.csv",
                    mapDict : dict = { 
                                        'LTR106' : 'LTR106_Mam' #nanopore ESCA
                                        },
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
    v35GeneMapDF = mapHUGO_2_ENSG(oldTxt2GeneFilePath)
    logger.info(f'v35GeneMapDF.shape : {v35GeneMapDF.shape }')
    selectGeneRows = v35GeneMapDF.loc[:, "HUGO"].isin(genes)
    logger.info(f'sum(selectGeneRows): {sum(selectGeneRows)}')
    v35RefBiomarkersDF = v35GeneMapDF.loc[selectGeneRows, :]
    logger.info(f'v35RefBiomarkersDF.shape : {v35RefBiomarkersDF.shape}')

    # elife transcripts counts where create using v39 repeat mask
    # load v39 mapping file and select genes of interest
    # i.e the v35 ENSG base == v39 ENSG base
    v39GeneMapDF= mapHUGO_2_ENSG(newTxt2GeneFilePath)
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

    retMapDF = mapDF.loc[:, ["HUGO_v35", "ENSG_v35", "ENSG_v39"]]
    if len(missingElifeGenes) > 0 :
        # we can not iterate and modify the same list
        tmpList = missingElifeGenes.copy()
        for oldName in tmpList:
            if oldName in mapDict:
                newName = mapDict[oldName]
                elifeGenes.append( newName )
                missingElifeGenes.remove( oldName )
                #add a row to mapDF
                retMapDF.loc[len(retMapDF.index)] = [oldName, newName, newName ] 
            else:
                logger.warning(f"missingElifeGenes: {missingInV35Set}")

    logger.info('END')
    return ( elifeGenes, missingElifeGenes, retMapDF )
