#
# selectSamples.py
# 
# Andrew E. Davidson
# aedavids@ucsc.edu
#

import logging
import pandas as pd

##############################################################################
class SelectSamples( object ):
    '''
    select and simplifies samples CIBERSORTx fractions results file

    public functions:
        select()
        simplify()

    ref: analysis/test/testSelectSamples.py
    '''

    logger = logging.getLogger(__name__)

    ###############################################################################
    def __init__( self,  
                    CIBERSORTxFractionsResultsFilePath : str,
                    separator : str = "\t"):
        '''
        TODO

        arguments:
            CIBERSORTxFractionsResultsFilePath : str
                example analysis/test/data/results.tsv

            separator :
                file type is expected to be tsv or csv.
                file column separator. Should be "," or "\t". 
                default = "\t"

        '''
        self.resultsFilePath = CIBERSORTxFractionsResultsFilePath
        self.sep = separator

        self.resultsDF = pd.read_csv(self.resultsFilePath, sep=self.sep, index_col="Mixture")

    ###############################################################################
    def select(self, mixtures : list[str]) -> pd.DataFrame :
        '''
        arguments:
            mixtures : list[str]
            list of sample ids found in the result file's 'Mixture' column 
        '''
        self.logger.info("BEGIN")
        retDF = self.resultsDF.loc[mixtures]

        self.logger.info("END")
        return retDF

    ###############################################################################
    def simplify(self, 
                df : pd.DataFrame,
                threshold : float = 0.1) -> pd.DataFrame :
        '''
        simplifies fractions by first setting the values less than threshold to zero.
        The selecting columns that do not sum to zero

        arguments:
            resultsDF:
                a pandas data frame in CIBERSORTx fractions results format

            threshold:
                default = 0.1
        
        '''
        self.logger.info("BEGIN")

        # we do not want to simpilfy any of the stats columns
        statCols = ['P-value', 'Correlation', 'RMSE']
        statsDF = df.loc[:, statCols]

        fractionsCols = ~df.columns.isin( statCols )
        self.logger.info(f'factionsCols: {fractionsCols}')
        fractionsDF = df.loc[:, fractionsCols]

        #
        # simplify
        #

        # where (): if cond is True, keep the original value. 
        # else replace with 0.0 
        thresholdDF = fractionsDF.where( fractionsDF >= threshold, 0.0)

        # do not select columns that sum to zero
        simplifedDF = thresholdDF.loc[ :, (thresholdDF.sum(axis=0) != 0) ]

        # add stats back
        byCols = 1
        retDF = pd.concat( [simplifedDF, statsDF], axis=byCols )

        self.logger.info("END")
        return retDF
