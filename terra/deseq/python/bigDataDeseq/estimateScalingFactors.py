'''
Created on Nov 8, 2021

@author: Andrew Davidson
aedavids@ucsc.edu

reference: prepareDataForDESeq2.ipynb
contains all the function need to calculate the estimated
'''

from   functools import reduce
import logging
from   operator import add
from   pyspark.sql.functions import col, exp, log, percentile_approx, round


###############################################################################
class EstimateScalingFactors( object ):
    '''
    DESeq is implemented in R and uses an R implementation detail to remove
    genes with zero in one or more samples. This will not work in spark. R defines
    log(0) = -inf. All arithmetic operations on -inf = -inf. DESeq calculates
    log of all values then calculates the row means. It then removes rows where
    row mean = -inf.

    In spark log(0) is defined to null. It seems like spark will treat null as
    zero. That is to say that if some genes have zero counts the rowSums will not
     be -inf. We will need to filter these genes out before calculating
    row means.

    a) calculate the log of all values
        i) logs are no easily swayed by outliers
        ii) we may have zero counts, log is undefined. These values will be null

    b) filter out genes with one or more nulls
       i) removes genes with zero in one or more samples.
        that are type specific. I.e. we want to focus on the house keeping genes.
        These are genes that are trascripted at similar levels regradless of
        tissue type

    c) calculate the 'geometic average' of the rows
        i) mean( rowsum )

    d) subtract the average log values from the log(counts)
        i) this is equal to log( numRead_x / average numRead_x)

    e) calculate the median of the ratio for each sample
        i) median is robust

    f) convert the medians back to linear scale
    '''

    logger = logging.getLogger( __name__ )

    ################################################################################
    def __init__( self, spark, txId2GeneIdFile, log4jLogger=None ):
        '''
        Constructor

        arguments:
            spark:
                an object return by pyspark.sql.SparkSession.builder

            # fileList:
            #     a file. each line is the path or url to salmon quant.sf file

            txId2GeneIdFile
                a csv file used to map transcript ids to gene ids. It has two columns with
                names 'txId'and 'geneId'

            log4jLogger
                The spark.sparkContext._jvm.org.apache.log4j
                default value is None ie. use python logging
        '''
        if log4jLogger is not None:
            self.logger = log4jLogger.LogManager.getLogger( __name__ )

        self.logger.warn( "__init__ BEGIN" )

        self.spark = spark
        # self.fileList = fileList
        # self.sampleNamesList = sampleNames
        self.txId2GeneIdFile = txId2GeneIdFile

        self.logger.warn( "__init__ END" )

    ################################################################################
    def run( self, rawCountsSparkDF, columnBatchSize, outFileGroupedByGene=None):
        '''
        Arguments:
            rawCountsSparkDF
                a dataframe with columns
                    The name column from the salmon quant.sf files

                    for each sample there is a column containing the NumReads
                    column of the salmon quant.sf file. The column name == the sample Name
                    
            columnBatchSize:
                an integer
                The GTEx training data set has 10409 numeric columns. This cause a
                java.lang.StackOverflowError because the DAG is to big. increasing spark driver
                memory does not help. The work around is sum the smaller batches of columns
                and cache the results of each batch
                
            outFileGroupedByGene:
                optional file path
                if defined groupedByGene dataframe will be saved. is a work around. 
                We are having trouble calculate the row sums needed to calculate the 
                estimated scaling factors. OOM exceptions. Using the grouped by gene counts matrix
                you can try having DESeq calculate the scaling factors with out having to 
                process all the salmon quant.sf file
                
        returns:
            (scalingFactorsDF, countDF)
                scalingFactorsDF:
                    a spark data frame with columns 'sampleName' and 'scalingFactor'

                countDF:
                    contains the integer counts of the transcripts grouped by geneId.
                    the first column name will be 'geneId'. The following column names will be the
                    sample names
        '''
        self.logger.warn( "run BEGIN" )
        
        self.logger.warn( "run rawCountsSparkDF numRows:{} numCols:{}"\
                         .format( rawCountsSparkDF.count(), len( rawCountsSparkDF.columns ) ) )
        

        # pass transients to enable unit testing
        rawCountsSparkDF.createOrReplaceTempView( "rawCounts" )

        countsSparkDF = self._groupByGeneAndSum( rawCountsSparkDF )
        retIntSparkDF = self._convertToLong( countsSparkDF )
        
        if outFileGroupedByGene:
            self.logger.warn("saving integer grouped by counts to :{}".format(outFileGroupedByGene))
            retIntSparkDF.coalesce(1).write.csv( outFileGroupedByGene, mode='overwrite', header=True)
            self.logger.warn("finished writing integer grouped by counts to :{}".format(outFileGroupedByGene))            
            
        countsSparkDF = None

        # 6.a)
        # skip first column, i.e. gene_id
        columnNames = retIntSparkDF.columns[1:]
        logCountsSparkDF = self._calculateLogs( retIntSparkDF, columnNames )

        # 6.b) filter out genes with one or more nulls
        #    i) removes genes with zero in one or more samples.
        #     that are type specific. I.e. we want to focus on the house keeping genes.
        #     These are genes that are trascripted at simpilar levels regradless of tissue type

        filteredDF = logCountsSparkDF.na.drop()
        filteredDF.checkpoint()
        logCountsSparkDF = None

        # 6.c) calculate the mean of the row sum
        # skip gene_id column
        columns = filteredDF.columns[1:]
        rowSumsDF = self.rowSums( filteredDF, columns, columnBatchSize )
        rowSumsDF.checkpoint()
        
        n = len( rowSumsDF.columns ) - 2  # do not count geneId or rowSum columns
        rowMeansDF = rowSumsDF.withColumn( "rowMean", ( rowSumsDF.rowSum / n ) )
        rowMeansDF.checkpoint()
        filteredDF = None

        # 6.d) subtract the avereage log values from the log(counts)
        #     i) this is equal to log( numRead_x / average numRead_x)

        # skip the first and last 2 columns, ie. geneId, rowSum, rowMean
        columnNames = rowMeansDF.columns[1:-2]
        ratioDF = self._subtractRowMeanFromLogCounts( rowMeansDF, columnNames )
        ratioDF.checkpoint()
        rowMeansDF = None

        # 6.e calculate the median of the ratio for each sample
        #     i) median is robust
        # skip geneId
        columnNames = ratioDF.columns[1:]
        logMedianDF = self.median( ratioDF, columnNames )
        logMedianDF.checkpoint()
        ratioDF = None

        newColNames = [self.getSampleNames( c ) for c in logMedianDF.columns]
        logScalingFactorsDF = logMedianDF.toDF( *newColNames )
        logScalingFactorsDF.checkpoint()

        # 6.f) convert the medians back to linear scale
        scalingFactorsDF = logScalingFactorsDF.select( *( exp( c )
                                     for c in logScalingFactorsDF.columns ) )

        # fix the column names change 'EXP(ctrl_1)' to 'ctrl_1'
        # transpose into a 2 columns, 'sampleName' and  'scalingFactor'
        retScalingFactorsDF = self._fixScalingFactors( scalingFactorsDF )

        # fix the column names change 'sum(kras)' to 'kras'
        retIntCountSparkDF = self._fixSumColNames( retIntSparkDF )

        self.logger.warn( "run END\n" )
        return ( retScalingFactorsDF, retIntCountSparkDF )

    ###############################################################################
    def _groupByGeneAndSum( self, rawCountSparkDF ):
        '''
        may change row ordering
        '''
        self.logger.warn( "_groupByGeneAndSum BEGIN" )

        rawCountSparkDF.createOrReplaceTempView( "rawCounts" )

        tx2geneSchema = " `txId` STRING, `geneId` STRING "
        txt2geneSparkDF = self.spark.read.load( self.txId2GeneIdFile, format="csv",
                                        schema=tx2geneSchema, header="false" )
        txt2geneSparkDF.createOrReplaceTempView( "txt2gene" )
        self.logger.warn( "_groupByGeneAndSum txt2geneSparkDF numRows:{} numCols:{}"\
                         .format( txt2geneSparkDF.count(), len( txt2geneSparkDF.columns ) ) )

        sqlStmt = 'select geneId, rc.* \n\
                        from \
                            rawCounts as rc, \n\
                            txt2gene  \n\
                        where \n\
                            rc.Name == txt2gene.txId'
        self.logger.debug( "sqlStmt:\n{}\n".format( sqlStmt ) )

        rawCountSparkDF = self.spark.sql( sqlStmt )
        rawCountSparkDF.createOrReplaceTempView( "rawCounts" )

        retSparkDF = rawCountSparkDF.groupBy( "geneId" ).sum()  # do not sort it does not matter .sort("geneId").show()

        self.logger.warn( "_groupByGeneAndSum END\n" )
        return retSparkDF

    ###############################################################################
    def _convertToLong( self, countsSparkDF ):
        self.logger.warn( "_convertToLong BEGIN" )

        # columns[1:], skip geneId col
        # https://mrpowers.medium.com/performing-operations-on-multiple-columns-in-a-pyspark-dataframe-36e97896c378
        # https://treyhunner.com/2018/10/asterisks-in-python-what-they-are-and-how-to-use-them/
        # countsSparkDF = countsSparkDF.select(*(col(c).cast("integer").alias(c) for c in countsSparkDF.columns[1:]))

        # round before casting to int to match R implementation

        roundSparkDF = countsSparkDF.select( col( 'geneId' ), \
                                              *( round( col( c ) ).alias( c ) for c in countsSparkDF.columns[1:] ) )

        self.logger.warn( "_convertToLong roundSparkDF numRows:{} numCols:{}"\
                         .format( roundSparkDF.count(), len( roundSparkDF.columns ) ) )

        intSparkDF = roundSparkDF.select( col( 'geneId' ), \
                                              *( col( c ).cast( "long" ).alias( c ) for c in countsSparkDF.columns[1:] ) )

        self.logger.warn( "_convertToLong intSparkDF numRows:{} numCols:{}"\
                         .format( intSparkDF.count(), len( intSparkDF.columns ) ) )

        self.logger.warn( "_convertToLong END\n" )
        return intSparkDF

    ###############################################################################
    def rowSums( self, countsSparkDF, columnNames, columnBatchSize ):
        '''
        The GTEx training data set has 10409 numeric columns. This cause a
        java.lang.StackOverflowError because the DAG is to big. increasing spark driver
        memory does not help. The work around is sum the smaller batches of columns
        and cache the results of each batch
        '''
        self.logger.warn("rowSums BEGIN")
        # tmpColumns = []
        totalColName = "rowSum"
        for i in range(0, len(columnNames), columnBatchSize) :
            tmpColName = "tmpSum" + str(i)
            # tmpColumns.append(tmpColName)
            batch = columnNames[i:i+columnBatchSize]
            countsSparkDF = self.rowSumsImpl(countsSparkDF, tmpColName, batch)
            
            if i == 0:
                countsSparkDF = countsSparkDF.withColumnRenamed(tmpColName, totalColName)
                             
            else:
                # calculate rolling total
                countsSparkDF = countsSparkDF.withColumn(totalColName, col(totalColName) + col(tmpColName))               
                # save space
                countsSparkDF = countsSparkDF.drop(tmpColName )                 
                
            # use an action to force execution
            numRows = countsSparkDF.count()
            self.logger.warn("rowSums:batch:{} numRows:{}".format(i, numRows))
            
            # check point will save the df data but not its linage
            countsSparkDF.checkpoint()                
        
            # self.logger.warn("AEDWIP remove show \n")
            # countsSparkDF.show()            
        
        
        # self.logger.warn("AEDWIP remove show final results")
        # countsSparkDF.show()
        
        self.logger.warn("rowSums END")
        return countsSparkDF       
    
    ###############################################################################
    def rowSumsImpl( self, countsSparkDF, newColName, columnNames ):
        '''
        calculates actual sum of columns
        
        arguments
            countSparkDF 
            
            newColumName: 
                results from column sum will be sorted here
            
            columnNames:
                list of columns to sum
                
        returns 
            amended countSparkDF 
        '''
        self.logger.warn( "rowSumsImpl BEGIN" )

        # https://stackoverflow.com/a/54283997/4586180
        retDF = countsSparkDF.na.fill( 0 ).withColumn( newColName , reduce( add, [col( x ) for x in columnNames] ) )

        # self.logger.warn( "rowSums retDF numRows:{} numCols:{}"\
        #                  .format( retDF.count(), len( retDF.columns ) ) )
        #
        # self.logger.warn("AEDWIP remove show")
        # retDF.show()

        self.logger.warn( "rowSumsImpl END\n" )
        return retDF

    ###############################################################################
    def _calculateLogs( self, countSparkDF, columns ):
        self.logger.warn( "_calculateLogs BEGIN" )

        colNameAndAlias = ( ( c, c.replace( 'sum', "log" ) )  for c in columns )
        ret = countSparkDF.select( col( 'geneId' ) ,
                                        *( log( c ).alias( a ) for c, a in colNameAndAlias ) )

        self.logger.warn( "_calculateLogs ret numRows:{} numCols:{}"\
                         .format( ret.count(), len( ret.columns ) ) )

        return ret
        self.logger.warn( "_calculateLogs END\n" )

    ###############################################################################
    def _subtractRowMeanFromLogCounts( self, sparkDF, columnName ):
        self.logger.warn( "_subtractRowMeanFromLogCounts BEGIN" )

        colNameAndAlias = ( ( c, c.replace( 'log(', '' ).replace( ')', '' ) )  for c in columnName )
        retDF = sparkDF.select( col( 'geneId' ),
                           *( ( col( c ) - col( 'rowMean' ) ).alias( a ) for c, a in colNameAndAlias ) )

        self.logger.warn( "_subtractRowMeanFromLogCounts retDF numRows:{} numCols:{}"\
                         .format( retDF.count(), len( retDF.columns ) ) )

        return retDF
        self.logger.warn( "_subtractRowMeanFromLogCounts END\n" )

    ###############################################################################
    def median( self, sparkDataFrame, columnNames ):
        '''
        https://spark.apache.org/docs/3.1.1/api/python/reference/api/pyspark.sql.functions.percentile_approx.html
        if you have an odd number of rows it does not calculate the value between the 2 middle values
        '''
        self.logger.warn( "median BEGIN" )
        retDF = sparkDataFrame.select( *( percentile_approx( c, 0.5, accuracy=1000000 )
                                         for c in columnNames ) )

        self.logger.warn( "median retDF numRows:{} numCols:{}"\
                         .format( retDF.count(), len( retDF.columns ) ) )

        self.logger.warn( "median END\n" )
        return retDF

    ###############################################################################
    def getSampleNames( self, colName ):
        # x = "percentile_approx(ctrl_1, 0.5, 1000000)"
        # name = x.split("(")[1].split(",")[0]
        ret = colName.split( "(" )[1].split( "," )[0]
    #     print(ret)
        return ret

    ###############################################################################
    def _fixSumColNames( self, df ):
        '''
        fix column names. the should be sample names
        # 'sum(ctrl_1)'' should be 'ctrl_1'
        '''
        self.logger.warn( "_fixSumColNames BEGIN" )

        # skip the the geneId column It is value like gene_2
        oldNames = df.columns[1:]

        sampleNames = ( ( c, c.replace( 'sum(', '' ).replace( ')', '' ) )  for c in oldNames )
        # print(*sampleColNames)

        retDF = df.select( col( 'geneId' ),
                               *( ( col( c ).alias( a ) for c, a in sampleNames ) ) )

        self.logger.warn( "_fixSumColNames END\n" )

        return( retDF )

###############################################################################
    def _fixScalingFactors( self, df ):
        '''
        df has 1 row, the scaling factors. the column names are 'EXP(ctrl_1)'
        change col names to 'ctrl_1', and transpose into 2 columns wiht names 'sampleName', and 'scalingFactor'
        '''
        self.logger.warn( "_fixScalingFactors BEGIN" )

        # # spark only method TODO test
        # https://newbedev.com/transpose-column-to-row-with-spark
        # from pyspark.sql import functions as func
        # #Use `create_map` to create the map of columns with constant
        # df = df.withColumn('mapCol', \
        #                     func.create_map(func.lit('col_1'),df.col_1,
        #                                     func.lit('col_2'),df.col_2,
        #                                     func.lit('col_3'),df.col_3
        #                                    )
        #                   )
        # #Use explode function to explode the map
        # res = df.select('*',func.explode(df.mapCol).alias('col_id','col_value')

        def getSampleNames( pdf ):
            # x = "EXP(ctrl_1)"
            columnNames = pdf.columns.to_list()
            ret = [ c.replace( "EXP(", "" ).replace( ")", "" ) for c in columnNames ]

            return ret

        # this should easily fit into driver memory
        # we want to save as 2 columns. It easier to transpose in pandas
        scalingFactorsPDF = df.toPandas()
        sampleNames = getSampleNames( scalingFactorsPDF )
    #     print("sampleNames:{}".format(sampleNames))
        scalingFactorsPDF = scalingFactorsPDF.transpose()

        # replace "EXP(sample name ) row names with sample name"
        scalingFactorsPDF = scalingFactorsPDF.set_axis( sampleNames, axis='index' )
        scalingFactorsPDF = scalingFactorsPDF.rename_axis( 'sampleName' ).reset_index()

        colNameDict = { 0:"scalingFactor"}
        colAxis = 1
        scalingFactorsPDF = scalingFactorsPDF.rename( colNameDict, axis=colAxis )

        retDF = self.spark.createDataFrame( scalingFactorsPDF )

        self.logger.warn    ( "_fixScalingFactors END\n" )
        return retDF
