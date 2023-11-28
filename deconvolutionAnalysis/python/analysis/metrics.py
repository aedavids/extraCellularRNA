#
# fractions.py
# functions for analysising the results from CIBERSORTx as if it was a multi-label classifier
#
# Andrew E. Davidson
# aedavids@ucsc.edu
#
# ref: extraCellularRNA/terra/jupyterNotebooks/cibersort/fractionsAsMulticlassClassification.ipynb
#

# FractionsCommandLine CommandLine display the doc string 
'''
    Treats the cibersort results as if they where the output layer of k-way (multiclass) classifier. The fractions can be thought of as a probablity distribution. This class choose use max to predict
    the class for each sample

    Calculates error metrics and statistics

    Do not confuse k-way classifiers with multilabel classifiers  
    
    See extraCellularRNA/deconvolutionAnalysis/python/analysis/test/testFractions.py
'''

import logging
import numpy as np
import os
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report

from analysis.metricsCLI import MetricsCommandLine

# global variables use by FractionsCommandLine
__all__ = []
__version__ = 0.1
__author__ = "Andrew Davidson aedavids@ucsc.edu"
__date__ = '2023-11-13'
__updated__ = '2023-11-13'

###############################################################################
class CibersortResultsAsKWayClassifier( object ):
    '''
    Treats the cibersort results as if they where the output layer of k-way (multiclass) classifier
    The fractions can be thought of as a probablity distribution. This class choose use max to predict
    the class for each sample

    Calculates error metrics and statistics

    Do not confuse k-way classifiers with multilabel classifiers

    Data overview:

    CIBERSORTX produce a results file. Each row contains the sampleId, the fractions, and stats
     Mixture, type1, type2, ... typeN, P-value, Correlation, RMSE

     for example
     Mixture
        GTEX-1117F-0226-SM-5GZZ7
        GTEX-1117F-0526-SM-5EGHJ
        GTEX-1117F-0726-SM-5GIEN

    ACC	Adipose_Subcutaneous	Adipose_Visceral_Omentum
        0.00000238818947786	0.32088468153413718	0.15663743975057184
        0.00002175019341662	0.29529977326057638	0.15036761048938591
        0.00015259963369839	0.05445776007930397	0.01055313063397182
        0.00000000000000000	0.29004114942525649	0.14226711911462919
        0.00021702689123660	0.00445898085087070	0.00000000000000000
        0.00012839158257415	0.00000000000000000	0.00000000000000000
        0.00017911589602714	0.06144797178936864	0.02300987639507408
        0.00000000000000000	0.00000000000000000	0.01792321073409516
        0.00073656971878856	0.20411159050680064	0.12949026476002565

    P-value	Correlation	RMSE
        0.00000000000000000	0.98542467644677711	0.92645530582020463
        0.00000000000000000	0.97979080799155893	0.93469518567034882
        0.00000000000000000	0.98490561378709340	0.44791559959147931
        0.00000000000000000	0.98846407133779690	0.90767508784039519
        0.00000000000000000	0.91705227835834124	0.94955257328411946
        0.00000000000000000	0.96999619583782004	0.71740096723763391
        0.00000000000000000	0.96897428384635731	0.80545057475847193
        0.00000000000000000	0.90218261892042595	1.23410866579020806
        0.00000000000000000	0.82951499979038523	0.95537827265096409
    
    '''
    logger = logging.getLogger(__name__)

    ###############################################################################
    def __init__( self,  
                 ciberSortResultsFilePath : str, expectedFractionsPath : str,
                 verbose : bool=False):
        '''

        arguments:
            ciberSortResultsFilePath : str
                see class documentation

            expectedFractionsPath : str
                path to file with same format as CIBESORTx results.txt file 
                with out the P-value, Correlation, RMSE columns.

                unlike the results.txt file the first column name is sampleId.


            verbose : bool
                default == false
                if true produce extra debug statements to stdout

        '''
        super().__init__()
        self.ciberSortResultsFilePath = ciberSortResultsFilePath
        self.expectedFractionsPath = expectedFractionsPath

        self.verbose = verbose

        # read data files
        rawCibersortFractionsDF, rawExpectedFractionsDF = self._load()
        self.rawCibersortFractionsDF = rawCibersortFractionsDF
        self.rawExpectedFractionsDF = rawExpectedFractionsDF

        self.fractionsWithPredictionsDF = self._getPredictionsDF( rawCibersortFractionsDF )
        self.labels = self._getAllClasses()
        ccc, ttt, ppp = self._getTrueAndPredictedTypes( self.fractionsWithPredictionsDF, 
                                                    rawExpectedFractionsDF )
        self.confusionDF = ccc
        self.trueLablesSeries = ttt
        self.predictedLabelsSeries = ppp

    ###############################################################################
    # def calculate(self, labels) -> tuple[np.array, str, dict]:
    def calculate(self, labels) -> tuple[np.array, pd.DataFrame]:
        '''
        Calculate metrics and statistics
        
        arguments:
            labels
                https://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html
                labels: List of labels to index the matrix. This may be used to reorder 
                or select a subset of labels. For evaluation metrics pass self.labels
            
                If None is given, those that appear at least once in y_true or y_pred are used in 
                sorted order. 

        returns 
            tuple(confusionMatrix, pandas.DataFrame)
            The dataframe contains stats and error metrics for each class

            see 
                https://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html#sklearn.metrics.confusion_matrix
                https://scikit-learn.org/stable/modules/generated/sklearn.metrics.classification_report.html#sklearn.metrics.classification_report
        '''
        # self.labels = sorted( set( self.confusionDF['category'].to_list() ) )
        # self.logger.info(f'unsortedlabels:\n{self.confusionDF["category"].to_list()}')
        # self.logger.info(f'sortedlabels:\n{self.labels}')

        # using the labels or not does not change the output of the confusion matrix?
        retConfusionMatrix = confusion_matrix(self.trueLablesSeries, 
                                                self.predictedLabelsSeries, 
                                                labels=labels)
        
        # https://scikit-learn.org/stable/modules/model_evaluation.html#classification-report
        # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.classification_report.html#sklearn.metrics.classification_report

       
        reportDict = classification_report(self.trueLablesSeries, 
                                        self.predictedLabelsSeries, 
                                        #target_names=labels, 
                                        labels=labels,
                                        zero_division=np.nan,
                                        output_dict=True)
        
        if not labels :
            # assume keys returned reportDict are something like
            # ['AAA', 'BBB', 'CCC', 'EEE', 'FFF', 'accuracy', 'macro avg', 'weighted avg']
            # we want all the 'class' keys
            classKeyList = list(reportDict.keys())
            classKeyList = classKeyList[:-3]
        else :
            classKeyList = labels
           
        sensitivitySpecificityDict = self._calculateCustomStats(retConfusionMatrix, classKeyList)

        # add sensitivity and specificity to final reportDict
        for className in classKeyList:
            classStatsDict = reportDict[className]
            extrStatsDict = sensitivitySpecificityDict[className]
            for key,value in extrStatsDict.items():
                classStatsDict[key] = value

        retDF = pd.DataFrame(reportDict).transpose()
        return (retConfusionMatrix, retDF)
    #
    # PRIVATE
    #

    ###############################################################################
    def _calculateCustomStats(self, confusionMatrix, labels) -> dict:
        '''
        ref: https://saturncloud.io/blog/how-to-report-sensitivity-and-specificity-in-multiclass-problems-using-scikitlearn/#sensitivity-and-specificity

        specificity = "true negative rate" = TN / (TN + FP)
        percision = "true positive rate" = TP / (TP + FP)
        sensitivity = recall = TP / (TP + FN)

        returns 
            a dictionary of dictionaries. for each class there is an inner dict. The keys are 
            specificity, sensitivity, tp, fn, fp, tn
            
            
            classification_report() returns recall and percision, allow redundent col to 
            make it easier to read

        '''
        numClasses = confusionMatrix.shape[0]

        retDict = {}
        for i in range(numClasses):
            className = labels[i]
            tp = confusionMatrix[i, i]
            fn = sum(confusionMatrix[i, :]) - tp
            fp = sum(confusionMatrix[:, i]) - tp
            tn = sum(sum(confusionMatrix)) - tp - fn - fp
            
            tpr = tp / (tp + fn)
            tnr = tn / (tn + fp)
            sensitivity = tp / (tp + fn)

            retDict[className] = {'specificity' : tnr, 
                                  'sensitivity' : sensitivity,
                                  'tp' : tp,
                                  'fn' : fn,
                                  'fp' : fp,
                                  'tn' : tn}

        return retDict
    
    ###############################################################################
    def _getAllClasses(self) -> list[str]:
        '''
        We need to get a list of all possible lables (ie classes) else it will be difficult
        to generate metrics and stats. self.confusionDF was contructed using an inner join
        and  may have dropped some class. 
        '''
        
        # the following commented code block does not work. 
        # we will not find classes for which we do not have any expected or predicted labels
        # expectedCategoriesNP = pd.unique(self.rawExpectedFractionsDF["category"])
        # predictedCategoriesNP = pd.unique(self.fractionsWithPredictionsDF['predictedCat'])
        # ret = np.sort( np.unique(np.concatenate( (expectedCategoriesNP, predictedCategoriesNP) )) )

        dropCols = ['participant_id', 'category', 'gender', 'age', 'dataSet']
        categorySet = set( self.rawExpectedFractionsDF.columns ) - set(dropCols)
        ret = sorted( list(categorySet) )

        self.logger.info(f"classes :{ret}")
        return ret

    ###############################################################################
    def _getPredictionsDF(self, rawCibersortFractionsDF : pd.DataFrame) -> pd.DataFrame :
        '''
        arguments
            rawCibersortFractionsDF : pd.DataFrame
            unfilterd data frame

        return : 
            -> pd.DataFrame
            copy of rawCibersortFractionsDF with a 'max' column. The value of the max column is the
            index into list of class types
        '''
        #
        # remove stat cols
        #
        statsCols = ['P-value', 'Correlation', 'RMSE' ]
        predictionsDF = rawCibersortFractionsDF.loc[:, ~rawCibersortFractionsDF.columns.isin( statsCols)].copy()

        #
        # treat predicition as if it was the output layer of a multiclassifier
        # use hard max to call the class
        #
        byRow = 1
        predictionsDF['max'] = predictionsDF.max(axis=byRow)

        #
        # find the name of the class associated with the max value
        # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.idxmax.html
        # Return index of first occurrence of maximum over requested axis
        #
        predictionsDF['predictedCat'] = predictionsDF.idxmax(axis=byRow)
        if self.verbose:
            print( f"predictionsDF.shape: {predictionsDF.shape}" )
        
        return predictionsDF
    ###############################################################################
    def _getTrueAndPredictedTypes( self,
                                    fractionsWithPredictionsDF : pd.DataFrame, 
                                    expectedFractionsDF : pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
        '''
        returns -> (pd.DataFrame, pd.Series, pd.Series)
            The DataFrame is the results of joining fractions and expected. We use join
            to insure that the data is order correctly.

            first series is the true category labels, the second is the predicted labels
        '''

        self.logger.debug(f'fractionsWithPredictionsDF\n{fractionsWithPredictionsDF}')
        self.logger.debug(f'expectedFractionsDF\n{expectedFractionsDF}')
        confusionDF = fractionsWithPredictionsDF.merge(right=expectedFractionsDF, 
                                                       how='inner', 
                                                       left_index=True, right_index=True,
                                                       suffixes=['pred', 'true'])
        # we have a lot of columns because both prediction and expected fractions have a col for each type
        if self.verbose :
            print( f"confusionDF.shape: {confusionDF.shape}" )

        return (confusionDF, confusionDF['category'], confusionDF['predictedCat'])


    
    ###############################################################################
    def _load(self) -> (pd.DataFrame, pd.DataFrame):
        '''
        reads 
        '''
        expectedFractionsDF  = pd.read_csv(self.expectedFractionsPath, sep='\t', index_col='sample_id')
        cibersortFractionsDF = pd.read_csv(self.ciberSortResultsFilePath, sep='\t', index_col='Mixture')

        if self.verbose:
            print( f"expectedFractionsDF.shape: {expectedFractionsDF.shape}" )
            print( f"cibersortFractionsDF.shape: {cibersortFractionsDF.shape}" )

        return (cibersortFractionsDF, expectedFractionsDF)

################################################################################
def main(inCommandLineArgsList=None):
    '''
    ref: pipeline/dataFactory/test/testInjection.py 

    ```
    $ python -m analysis.metrics -e analysis/test/data/expectedFractions.tsv \
                                 -f analysis/test/data/results.tsv \
                                 -o ./tmp
    ```

    see unit test unit test 
    extraCellularRNA/deconvolutionAnalysis/python/analysis/test/testFractions.py
    
    '''
    # we only configure logging in main module
    loglevel = "WARN"
    #loglevel = "INFO"
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
    cli = MetricsCommandLine( 
                    version=__version__ , 
                    author=__author__ ,
                    date=__date__, 
                    update=__updated__)

    if inCommandLineArgsList is None:
        cli.parse()
    else:
        cli.parse( inCommandLineArgsList )

    logger.warning(f'command line arguments : {cli.args}')


    fractionsPath         = cli.args.fractionsPath
    expectedFractionsPath = cli.args.expectedFractionsPath
    outDir                = cli.args.outDir

    os.makedirs(outDir, exist_ok=True)

    evaluator = CibersortResultsAsKWayClassifier(fractionsPath, expectedFractionsPath, verbose=True)

    #
    # find all possible labels. ie categories,tissue type, cell types, ... 
    # these are the we constructed the signature matrix from
    # if lables=None we only get back metric for classes that have at least one y_pred 
    # y_true. We want to make sure we can identify classes for which we had not samples
    #

    confusionMatrixNP, metricsDF = evaluator.calculate(labels=evaluator.labels)
    metricsDF.index.name = "id"

    confusionMatrixPath = os.path.join(outDir, "confusionMatrix.csv")
    # np.savetxt(confusionMatrixPath, confusionMatrixNP, delimiter=",", fmt="%s")
    confusionMatrixDF = pd.DataFrame(confusionMatrixNP, 
                                     index=evaluator.labels,
                                     columns=evaluator.labels,
                                     )
    confusionMatrixDF.index.name="true\predicted"
    confusionMatrixDF.to_csv(confusionMatrixPath)
    logger.warning(f'saved confussion matrix to : {confusionMatrixPath}')

    metricsPath = os.path.join(outDir, "metrics.csv" )
    logger.debug(f'AEDWIP metricsDF\n {metricsDF}')
    #metricsDF.round(decimals=3).to_csv(metricsPath)
    metricsDF.to_csv(metricsPath)
    logger.warning(f'saved metrics to : {metricsPath}')



################################################################################
if __name__ == '__main__':
    main()
