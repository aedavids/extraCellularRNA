#
# createCibersortSignatureMatrix.py
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref:
# extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
# extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/cibersortMixtureMatrixFactory.py

# createCiberSortInputMatrices cli the doc string 
'''
create signatureMatrix.tsv file for cbersortx from count data created by 
"Create", formally known as 'Complete Seq'
    
ref:
    extraCellularRNA/terra/jupyterNotebooks/cibersort/createCiberSortGeneSignatureMatrix.ipynb
    extraCellularRNA/deconvolutionAnalysis/python/pipeline/dataFactory/cibersortMixtureMatrixFactory.py
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from kimLabUtils.baseCommandLine import BBaseCommandLine
import logging
logger = logging.getLogger(__file__)
import os
import pandas as pd
import sys

__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2024-12-23'
__updated__ = '2024-12-23'

###############################################################################
class CreateSignatureMatrixCommandLine( BBaseCommandLine ):
    '''
    Handle the command line, usage and help requests.
    '''
    ###############################################################################
    def __init__( self, version, author, date, update ):
        '''
        Implement a parser to interpret the command line argv string using argparse.

        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''
        super().__init__(version, author, date, update )
        # self.author = author
        # self.program_version = version
        # self.program_build_date = str( update )
        # self.date = date
        #
        # self.program_version_message = '%%(prog)s %s (%s)' % ( self.program_version, self.program_build_date )
        # self.program_shortdesc = __import__( '__main__' ).__doc__.split( "\n" )[1]

    ###############################################################################
    def _build( self ):
        self.parser = ArgumentParser( description=self._getLicence(), formatter_class=RawDescriptionHelpFormatter )
    
        #
        # optional arguments
        #

        # self.requiredArg.add_argument( '-g', '--genesOfInterest', required=True, nargs='+', metavar="",
        #                           action='store', 
        #                           help="list of genes of interest"
        # )
        # --categoriesOfInterest
        self.parser.add_argument( '-c', '--categoriesOfInterest', default=None, required=False, action='store', nargs='+', metavar="",
                                  help="list of categories to include in the signature matrix. If None include all categories."
                                        + " Example of use: You want to a limit of detection or dilution fractions series."
                                        + " Your categories are control, undiluted, 1:10, 1:100, 1:1000. Your categoriesOfInterest"
                                        + "would be ['control', 'undiluted']"
        )
       
        self.parser.add_argument( '-u', '--useMedian', default=False, action='store_true', 
                      help="use median to calculate the signature gene mean value for each category. Default is to use the mean"
        )

        # self.parser.add_argument( '-o', '--outDir', default=".", metavar="",
        #                                               action='store', 
        #                                               help="locaiton to write output files"
        # )

    
        #
        # group required arguments. This will create a better help message
        # make sure to set required=True
        #
        self.requiredArg = self.parser.add_argument_group( 'required arguments' )

        self.requiredArg.add_argument( '-o', '--outDir', default=".", metavar="",
                                                      action='store', 
                                                      help="locaiton to write output files mixture.tsv and signatureGenes.tsv"
        )        

        self.requiredArg.add_argument( '-n', '--normalizedCountFilePath', required=True, default=None, metavar="",
                                                      action='store', 
                                                      help="path to a csv file with normalized gene counts. "
                                                        + "ex. '/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv'"
                                                        + "1st 3 column are gene_id, gene_name, gene_biotype there is a column for each sample"
                                    )
        
        self.requiredArg.add_argument( '-m', '--metaDataFilePath', required=True, default=None, metavar="",
                                                      action='store', 
                                                      help="path to a csv file containing sample meta data in DESeq format"
                                                      + "ex. /private/groups/kimlab/data/tempus/illumina/20241107/raw/metaDataWithHeader.csv"
                                                      + "this file should not have a header, the first column should be sampleId, the second the sample type"
        )

        self.requiredArg.add_argument( '-g', '--genesOfInterest', required=True, nargs='+', metavar="",
                                  action='store', 
                                  help="list of genes of interest"
        )

################################################################################
def createSignatureMatrix(
        genesOfInterest : list[str],
        normalizedCountDF : pd.DataFrame,
        metaDF : pd.DataFrame,
        useMedian : bool = False,
        categoriesOfInterest : list[str] = None
    ) :
    '''
    create a signature matrix for CiberSort
    
    arguments:
        genesOfInterest : list[str]

        normalizedCountDF : pd.DataFrame
            The DESeq2 normalized count data. Each row is a gene.
            the index column should be named 'gene_id'

            sample data:
                gene_id, gene_name, gene_biotype, SLDK3_T1_100_S1_L007
                (A)n, (A)n, Microsatellite, 0.5450089567650692
                (AAA)n, (AAA)n, Microsatellite, 0
                (AAAAAAC)n, (AAAAAAC)n, Microsatellite, 0
                (AAAAAAG)n, (AAAAAAG)n, Microsatellite, 0

        metaDF : pd.DataFrame,
            the DESeq2 colData. This data describes each sample. 
            E.G. sample_id category, ...
            The index column should be note be named 'sample_id' or 'category'

            ./illumina/20241107/raw/metadata.csv

            sample data:
                sample_id              category
                SLDK3_T1_100_S1_L007   A
                SLDK3_T1_100K_S2_L007  B
                SLDK3_T1_10K_S3_L007   A
                SLDK3_T1_1K_S4_L007    B

        useMedian : bool
            if True use the median to calculate the signature gene mean value for each category
            if False use the mean to calculate the signature gene mean value for each category

        categoriesOfInterest : list[str]
            list of categories to include in the signature matrix. If None include all categories.
            Example of use: You want to a limit of detection or dilution fractions series.
            Your categories are control, undiluted, 1:10, 1:100, 1:1000. Your categoriesOfInterest
            would be ['control', 'undiluted']
            
        
    returns:
        pd.DataFrame with index name 'gene_id'. the column names are the list of categories

    ref:
        extraCellularRNA/deconvolutionAnalysis/python/tempus/test/testCreateCiberSortInputMatrices.py
    '''
    logger.info("BEGIN")
    logger

    # sort the gene list and remove any duplicates
    logger.info(f'normalizedCountDF.shape : {normalizedCountDF.shape}')   
    logger.info(f'metaDF.shape : {metaDF.shape}')   


    transposedDF = normalizedCountDF.transpose(copy=True)
    logger.info(f'transposedDF.shape : {transposedDF.shape}')   
    logger.debug(f'transposedDF\n{transposedDF}')   

    # join the metaDF, we need the 'category' col so we can
    # calculate the  signature gene mean value for each category 
    joinDF =  pd.merge(left=transposedDF, 
                        right=metaDF.loc[:,["sample_id", "category"]], 
                        how='inner', 
                        left_index=True, 
                        right_on="sample_id")      
    logger.info(f'joinDF.shape : {joinDF.shape}')
    logger.info(f'joinDF.head():\n{joinDF.head()}')

    if categoriesOfInterest is not None:
        selectRows = joinDF.loc[:, "category"].isin(set(categoriesOfInterest))
        joinDF = joinDF.loc[selectRows, :]

    # calculate the expected values for each category  
    sortedGeneList = sorted( set(genesOfInterest) ) 
    genesDF = joinDF.loc[ :, sortedGeneList + ["category"] ] 
    logger.info(f'genesDF.shape : {genesDF.shape}')

    if useMedian :
        # weird duplicated log so I can set debugger break points
        logger.info(f'useMedian : {useMedian} calling median()')
        signatureDF = genesDF.groupby("category").median()
    else:
        logger.info(f'useMedian : {useMedian} calling mean()')
        signatureDF = genesDF.groupby("category").mean()
    
    # convert to cibersort expected upload format
    ciberSortSignatueDF = signatureDF.transpose(copy=True)
    logger.info(f'ciberSortSignatueDF.shape : {ciberSortSignatueDF.shape}')
    
    ciberSortSignatueDF.index.name = "gene_id"
    
    # weird. cciberSortSignatueDF.columns.name = 'category'. This name is not
    # saved by pd.to_csv(). Set to none to make it easier to write unit test
    ciberSortSignatueDF.columns.name = None

    logger.info(f'END')
    return ciberSortSignatueDF


################################################################################
def main(inCommandLineArgsList=None):
    '''
    TODO
    '''
    # we only configure logging in main module
    # loglevel = p.getProperty("LOG_LEVEL")
    #loglevel = "INFO"
    loglevel = "WARN"
    # logFMT = p.getProperty("LOG_FMT")
    # logFMT = "%(asctime)s %(levelname)s [thr:%(threadName)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logFMT = "%(asctime)s %(levelname)s %(name)s %(funcName)s() line:%(lineno)s] [%(message)s]"
    logging.basicConfig(format=logFMT, level=loglevel)    

    logger = logging.getLogger(__file__)

    logger.warning("BEGIN")

    #
    # always log run time env to make debug easier
    #
    ORIG_PYTHONPATH = os.environ['PYTHONPATH']
    logger.warning(f'PYTHONPATH : {ORIG_PYTHONPATH}')
    logger.warning(f'FILE: {__file__}')
    logger.warning(f'PWD: {os.getcwd()}')

    #
    # parse the command line
    #
    cli = CreateSignatureMatrixCommandLine( 
        version=__version__,
        author=__author__, 
        date=__date__,
        update=__updated__ )

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')

    outDir = cli.args.outDir
    normalizedCountFilePath = cli.args.normalizedCountFilePath
    metaDataFilePath = cli.args.metaDataFilePath
    genesOfInterest = cli.args.genesOfInterest
    categoriesOfInterest = cli.args.categoriesOfInterest

    if cli.args.useMedian:
        useMedian = True
    else :
        useMedian = False
    logger.info(f' useMedian : {useMedian}')

    # make sure output director exist
    os.makedirs(outDir, exist_ok=True)

    #
    # load the normalized counts and get rid of any extra columns
    #
    normalizedCountDF = pd.read_csv(normalizedCountFilePath)

    normalizedCountDF.set_index("gene_id", inplace=True)
    
    if "gene_name" in normalizedCountDF.columns:
        normalizedCountDF.drop(columns=["gene_name"], inplace=True)

    if "gene_biotype" in normalizedCountDF.columns:
        normalizedCountDF.drop(columns=["gene_biotype"], inplace=True)

    logger.info(f'normalizedCountDF.shape : {normalizedCountDF.shape}')

    #
    # load the meta data
    #
    metaDataDF = pd.read_csv(metaDataFilePath)
    logger.info(f'metaDataDF.shape : {metaDataDF.shape}')
    logger.info(f'metaDataDF :\n {metaDataDF}')

    #
    # create the cibersort signature matrix
    #
    retDF = createSignatureMatrix(
            genesOfInterest, 
            normalizedCountDF,
            metaDataDF,
            useMedian,
            categoriesOfInterest
        )

    logger.info(f'retDF.shape():\n{retDF.shape}')

    # save output file
    os.makedirs(outDir, exist_ok=True)
    outPath = f'{outDir}/signatureMatrix.tsv'
    
    retDF.to_csv(outPath, index=True, sep='\t')

    print(f'saved signature matrix to {outPath}')

    logger.warning("END")
    sys.exit(0)

################################################################################
if __name__ == '__main__':
    # debugCommandLineArgsList=[
    #     #"--help",
    #     #"-useMedian",
    #     "--outDir", "./tmp",
    #     "--normalizedCountFilePath", "/private/groups/kimlab/data/tempus/illumina/20241107/create/annotated_norm_counts.csv",
    #     "--metaDataFilePath", "/private/groups/kimlab/data/tempus/illumina/20241107/raw/metaDataWithHeader.csv",
    #     "--genesOfInterest", "X7D_LINE", "Zaphod", 
    #     "--categoriesOfInterest", "Control", "UD"
    # ]    
    # main(debugCommandLineArgsList)

    main()
