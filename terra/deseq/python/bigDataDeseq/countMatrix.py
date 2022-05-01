'''
Created on Jan 18, 2022

@author: andrewdavidson
'''

import logging
# from   pyspark.sql.functions import *
# from   pyspark.sql.window import *

###############################################################################
class CountMatrix(object):
    '''
    classdocs
    '''

    logger = logging.getLogger( __name__ )

    ###############################################################################
    def __init__(self, spark, log4jLogger=None):
        '''
        arguments
            spark:
                an object return by pyspark.sql.SparkSession.builder
        
           log4jLogger
                The spark.sparkContext._jvm.org.apache.log4j
                default value is None ie. use python logging
        '''
        if log4jLogger is not None:
            self.logger = log4jLogger.LogManager.getLogger( __name__ )

        self.logger.info( "BEGIN" )
        self.spark = spark
        
        self.logger.info( "END" )

    ###############################################################################    
    def _addRowIdx(self, df, colName):
        '''
        # https://kb.databricks.com/sql/gen-unique-increasing-values.html
        '''
        # window = Window.orderBy(col('monotonically_increasing_id'))
        # df_with_consecutive_increasing_id = df_with_increasing_id.withColumn('increasing_id', row_number().over(window))
        # df_with_consecutive_increasing_id.show()
        #
        # retDF = df.withColumn('index_column_name', row_number().over(Window.orderBy(monotonically_increasing_id())) - 1) 
        
    
        df = df.withColumn("mono_increase_idx", monotonically_increasing_id() )
        window = Window.orderBy( col("mono_increase_idx") )
        retDF = df.withColumn( colName, row_number().over(window) - 1)
        
        return retDF

    ###############################################################################    
    def loadSalmonReadsTable(self, fileList, sampleNamesList):
        '''
        for each file in file list, selects a spark data frame with the 'name column'
        and the numReads column renamed as the sample name
        
        The data frames are then joined together where quant file names are equal 
        into a single count dataframe. 
        
        WARNING: this does not work well with large data set stored in txt,CSV,TSV format. 
        Name col with a range value. Name col + numReads is 491 m. on disk
        just numReads is 21 m on disk SEE: loadSalmonReadsTableWithRowId(). 
        Try using Data Lake/Parquet file format
        
        arguments
        
                fileList:
                    a file. each line is the path or url to salmon quant.sf file

                sampleNameslist
                    must be in same order as fileList
        '''
        self.logger.info( "BEGIN" )
        retNumReadsDF = None
        quantSchema = "`Name` STRING, `Length` INT, `EffectiveLength` DOUBLE, `TPM` DOUBLE, `NumReads` DOUBLE "
        for i in range( len(fileList) ):
            #
            # get NumReads from next salmon quant file
            #
            quantFile = fileList[i]
            sampleDF = self.spark.read.load( quantFile, format="csv", sep="\t",
                                     schema=quantSchema, header="true" ) 
                                    # did not fix bug .repartition(50)

            sampleName = sampleNamesList[i]
            sampleDF = sampleDF.select( ["Name", "NumReads"] )\
                            .withColumnRenamed( "NumReads", sampleName )
            
            sampleDF.createOrReplaceTempView( "sample" )
                        
            self.logger.debug("AEDWIP i:{} sampleName:{} sampleDF.num rows:{} num cols:{} num parts:{}"
                             .format(i, sampleName, sampleDF.count(), len(sampleDF.columns), sampleDF.rdd.getNumPartitions()))
            
            #
            # append NumReads to table of reads
            #
            
            # the sample name must be quoted else column names with a '-'
            # like GTEX-1117F-0426-SM-5EGHI will generate an error
            # spark think the '-' is an expression. '_' is also
            # a special char for the sql like operator
            # https://stackoverflow.com/a/63899306/4586180
            sqlStmt = '\t\t\t\t\t\tselect rc.*, `{}` \n\
                            from \n\
                               retNumReadsDF as rc, \n\
                               sample  \n\
                            where \n\
                                rc.Name == sample.Name \n'.format( sampleName )

            # funky pythonic way to write goal print a nice string
            # sqlStmt =(f"select rc.*, {sampleName} ",
            #         "from rawCounts as rc, sample",
            #         " where ",
            #         "rc.Name == sample.Name")
            # self.logger.warn(
            # sqlStmt = sqlFmt.format(sampleName)

            self.logger.debug( "sqlStmt:\n{}\n".format( sqlStmt ) )
            if i == 0 :
                retNumReadsDF = sampleDF
            else :
                retNumReadsDF = self.spark.sql( sqlStmt )
                
            retNumReadsDF.createOrReplaceTempView( "retNumReadsDF" )
            
            #
            # debug. seems like we do not make progress when we run on GTEx training
            # nothing happens, logs do not change, cluster metrics drop suggesting no work
            # is being done
            # add an action to try and debug
            # this should not change the physical plan. I.e. we still have the same number of shuffles
            # which results in the same number of stage. We are just not building up a plan with thousands
            # of stages. 
            #
            self.logger.warn("AEDWIP i:{} retNumReadsDF.num rows:{} num cols:{} num parts:{}"
                             .format(i, retNumReadsDF.count(), len(retNumReadsDF.columns), retNumReadsDF.rdd.getNumPartitions()) )

            #
            # TODO AEDWIP spark analyze chapter 18 debugging joins
            
            # execution plan should be the same for each join
            #rawCountsSDF.explain()

        self.logger.info( "END\n" )
        return retNumReadsDF  
        

        
    ###############################################################################    
    def loadSalmonReadsTableWithRowId(self, fileList, sampleNameList):
        '''
        22/01/19 14:02:19 WARN WindowExec: No Partition Defined for Window operation! Moving all data to a single partition, this can cause serious performance degradation.

        for each file in file list, selects a spark data frame with the 'name column'
        and the numReads column renamed as the sample name
        
        The data frames are then joined together where quant file names are equal 
        into a single count dataframe. 
        
        arguments
        
                fileList:
                    a file. each line is the path or url to salmon quant.sf file

                sampleNamesList:
                    a list must be in same order as fileList        
        '''
        self.logger.info( "BEGIN" )
        retNumReadsDF = None
        quantSchema = "`Name` STRING, `Length` INT, `EffectiveLength` DOUBLE, `TPM` DOUBLE, `NumReads` DOUBLE "
        namesDF = None
        # numRows = 0
        
        for i in range( len(fileList) ):
            #
            # get NumReads from next salmon quant file
            # create an id we can join on. Names is 470 mb on disk
            # numReads is about is 21 m on disk
            sampleName = sampleNameList[i]
            quantFile = fileList[i]
            sampleDF = self.spark.read.load( quantFile, format="csv", sep="\t",
                                     schema=quantSchema, header="true" ) 
            
            idxColName = "rowIdx"
            if i == 0:
                # get a copy of the salmon quant.sf file's Names col
                namesDF = sampleDF.select( ["Name"] )  
                namesDF = self._addRowIdx(namesDF, idxColName)
                
            sampleDF = sampleDF.select( ["NumReads"] )\
                        .withColumnRenamed( "NumReads", sampleName )
            sampleDF = self._addRowIdx(sampleDF, idxColName)            

            
            self.logger.debug('AEDWIP REMOVE SHOW()')
            # sampleDF.show()
            
            sampleDF.createOrReplaceTempView( "sample" )
                        
            self.logger.debug("AEDWIP i:{} sampleName:{} sampleDF.num rows:{} num cols:{} num parts:{}"
                             .format(i, sampleName, sampleDF.count(), len(sampleDF.columns), sampleDF.rdd.getNumPartitions()))
            
            #
            # append NumReads to table of reads
            #
            
            # the sample name must be quoted else column names with a '-'
            # like GTEX-1117F-0426-SM-5EGHI will generate an error
            # spark think the '-' is an expression. '_' is also
            # a special char for the sql like operator
            # https://stackoverflow.com/a/63899306/4586180
            sqlStmt = '\t\t\t\t\t\tselect rc.*, `{}` \n\
                            from \n\
                               retNumReadsDF as rc, \n\
                               sample  \n\
                            where \n\
                                rc.{} == sample.{} \n'.format( sampleName, idxColName, idxColName )

            # funky pythonic way to write goal print a nice string
            # sqlStmt =(f"select rc.*, {sampleName} ",
            #         "from rawCounts as rc, sample",
            #         " where ",
            #         "rc.Name == sample.Name")
            # self.logger.warn(
            # sqlStmt = sqlFmt.format(sampleName)

            self.logger.debug( "sqlStmt:\n{}\n".format( sqlStmt ) )
            if i == 0 :
                retNumReadsDF = sampleDF
            else :
                retNumReadsDF = self.spark.sql( sqlStmt )
                
            retNumReadsDF.createOrReplaceTempView( "retNumReadsDF" )
            
            #
            # debug. seems like we do not make progress when we run on GTEx training
            # nothing happens, logs do not change, cluster metrics drop suggesting no work
            # is being done
            # add an action to try and debug
            # this should not change the physical plan. I.e. we still have the same number of shuffles
            # which results in the same number of stage. We are just not building up a plan with thousands
            # of stages. 
            #
            self.logger.debug("AEDWIP i:{} retNumReadsDF.num rows:{} num cols:{} num parts:{}"
                             .format(i, retNumReadsDF.count(), len(retNumReadsDF.columns), retNumReadsDF.rdd.getNumPartitions()) )

            #
            # TODO AEDWIP spark analyze chapter 18 debugging joins
            
            # execution plan should be the same for each join
            #rawCountsSDF.explain()
            
        #
        # add the 'Names' column back in and remove the 'tid' column and 
        #
        retNumReadsDF.createOrReplaceTempView( "retNumReadsDF" )
        namesDF.createOrReplaceTempView( "names" )
        sqlStmt = 'select names.Name, rc.* \n\
                            from \n\
                               retNumReadsDF as rc, \n\
                               names  \n\
                            where \n\
                                rc.{} == names.{} \n'.format(idxColName, idxColName)

        retNumReadsDF = self.spark.sql( sqlStmt )
        retCols = [c for c in retNumReadsDF.columns if c != idxColName and c != "mono_increase_idx" ]
        retNumReadsDF = retNumReadsDF.select( retCols )
        
        self.logger.info( "END\n" )
        return retNumReadsDF
    
    ###############################################################################
    def loadSalmonQuantFiles( self ):
        '''
        return:
            a list of spark data frames, one for each salmon quant.sf file
        '''
        self.logger.info( "BEGIN" )

        # pre allocate slots to store data frames in
        quantSDFs = [None] * len( self.fileList )

        quantSchema = "`Name` STRING, `Length` INT, `EffectiveLength` DOUBLE, `TPM` DOUBLE, `NumReads` DOUBLE "

        for i in range( len( quantSDFs ) ):
            quantFile = self.fileList[i]
            self.logger.debug( "quantFile:{}".format( quantFile ) )

            df = self.spark.read.load( quantFile, format="csv", sep="\t",
                                     schema=quantSchema, header="true" )
            quantSDFs[i] = df

        self.logger.info( "END\n" )
        return quantSDFs
    