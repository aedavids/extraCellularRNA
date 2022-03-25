'''
Created on Dec 22, 2021

@author: andrewdavidson
'''
import logging

class ParseSalmonReads(object):
    '''
    classdocs
    '''
    logger = logging.getLogger( __name__ )

    ################################################################################
    def __init__(self, spark, fileList, sampleNames, log4jLogger=None ):
        '''
        Constructor

        arguments:
            spark:
                an object return by pyspark.sql.SparkSession.builder

            fileList:
                a file. each line is the path or url to salmon quant.sf file
                
            sampleNames:
                a parallel list of sample names. i.e. sampleNames[0] is the 
                sample name for fileList[0]
                
            log4jLogger
                The spark.sparkContext._jvm.org.apache.log4j
                default value is None ie. use python logging
        '''
        if log4jLogger is not None:
            self.logger = log4jLogger.LogManager.getLogger(__name__)
            self.logger.warn("AEDWIP using log4j")
            
        self.logger.info( "BEGIN" )
        
        self.spark = spark
        self.fileList = fileList
        self.sampleNamesList = sampleNames
        
        self.logger.info( "END\n" )        
    
    ###############################################################################
    def selectCounts_depercated( self ):
        '''
        return:
            a list of spark data frames, one for each salmon quant.sf file
            each data frame as 2 columns 'name' and 'numReads'
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
            
            sampleName = self.sampleNamesList[i]
            df = df.select( ["Name", "NumReads"] )\
                    .withColumnRenamed( "NumReads",  sampleName)
            
            quantSDFs[i] = df

        self.logger.info( "END\n" )
        return quantSDFs


    ###############################################################################
    def selectCounts( self, outputDir, selectOnlyNumReads=False ):
        '''
        return:None
        
        TODO: a clean design make writing unit test easier is do not use for loop. or data members
            1. pass in df
            2. load
            3. select
            4 return df
            
            have main cli do write
           
        '''
        self.logger.info( "BEGIN" )

        quantSchema = "`Name` STRING, `Length` INT, `EffectiveLength` DOUBLE, `TPM` DOUBLE, `NumReads` DOUBLE "
        for i in range( len( self.fileList ) ):
            quantFile = self.fileList[i]
            self.logger.warn( "i:{} sampleName:{}".format(i, self.sampleNamesList[i] ) )
            self.logger.warn( "i:{} quantFile:{}".format(i, quantFile ) )

            df = self.spark.read.load( quantFile, format="csv", sep="\t",
                                     schema=quantSchema, header="true" )
            
            sampleName = self.sampleNamesList[i]
            
            if selectOnlyNumReads:
                colList = ["NumReads"] 
            else:
                colList = ["Name", "NumReads"] 
                
            df = df.select( colList )\
                    .withColumnRenamed( "NumReads",  sampleName)
            
            outfile = outputDir + "/" + sampleName 
            msg = "AEDWIP outfile:{}".format( outfile )
            self.logger.warn( msg )
            df.write.csv( outfile, mode='overwrite', header=True )
            df.unpersist()
        self.logger.info( "END\n" )

    ###############################################################################
    def split_deprecated(self, df, numberOfSplits):
        '''
        df is ordered by 'Name' then divides dataframe into numberOfSplits by row.
        the number of rows will be df.count() // numberOfSplits
        
        Note if count is no divisible by numberOfSplits the remain rows will be
        return in an extra dataframe
        '''
        self.logger.info( "BEGIN" )
        # Calculate count of each dataframe rows
        count = df.count()
        numRows = count // numberOfSplits
        numRows = int(numRows)
        remainingRows = count % numberOfSplits
        
        self.logger.warn("count:{} numberOfSplits:{} numRows:{} remainingRows:{}"
                         .format(count, numberOfSplits, numRows, remainingRows))
        
        # pre allocate results buffer
        if remainingRows > 0:
            retList = [None] * (numberOfSplits + 1)
        else :
            retList = [None] * numberOfSplits            


        # Create a copy of original dataframe
        copyDF = df.orderBy("Name")
        # copyDF.show()
        
        i = 0
        while i < numberOfSplits:
            self.logger.warn("i:{}".format(i))
            # Get the top `numRows` number of rows
            # note take() is an action
            # limit() is a transformation
            topDF = copyDF.limit( numRows )
            
            # Truncate the `copy_df` to remove
            # the contents fetched for `temp_df`
            # original quant.sf files are sorted by name however
            # we must use order by, else the row names between
            # GTEx sample will not be the same
            # we can not simply sort or orderBy once. we have to 
            # do this on every iteration
            copyDF = copyDF.subtract(topDF).orderBy( "Name" )

            retList[i] = topDF
            
            # Increment the split number
            i += 1
            
        if remainingRows > 0 :
            self.logger.info("AEDWIP writing last i:{} len(retList):{}".format(i, len(retList)))
            retList[i] = copyDF      
            #copyDF.show()
            #retList[i].show()
            
            
        self.logger.warn("AEDWIP debug split() retlist[0].count(): {}".format(retList[0].count()))
        
        self.logger.info( "END\n" )
        return retList

    ###############################################################################
    def writeCSV(self, df, logTag, filePath):
        msg = "{} filePath:{}".format( logTag, filePath )
        self.logger.warn( msg )
        df.write.csv( filePath, mode='overwrite', header=True )
        
    ###############################################################################
    def split(self, df, numberOfSplits, sampleName, outputDir, logTag):
        '''
        df is ordered by 'Name' then divides dataframe into numberOfSplits by row.
        the number of rows will be df.count() // numberOfSplits
        
        Note if count is no divisible by numberOfSplits the remain rows will be
        return None
        
        This version calls writeCSV after each split to add in debugging
        
        TODO: orderBy is expensive would be better to us windowing
        '''
        
        
#         aedwip see
#         https://lists.apache.org/thread/c31kxt8yvgt4m49tt35p5bcyy60bxhoq
#         lease try repartitionbyrange,
# dpark 3 has adaptive query execution with configurations to handle skew as
# well.
#
# Regards,
# Gourav        
        self.logger.info( "BEGIN" )
        # Calculate count of each dataframe rows
        count = df.count()
        numRows = count // numberOfSplits
        numRows = int(numRows)
        remainingRows = count % numberOfSplits
        
        self.logger.warn("count:{} numberOfSplits:{} numRows:{} remainingRows:{}"
                         .format(count, numberOfSplits, numRows, remainingRows))         

        # Create a copy of original dataframe
        copyDF = df.orderBy("Name")
        # copyDF.show()
        
        i = 0
        while i < numberOfSplits:
            self.logger.warn("i:{}".format(i))
            # Get the top `numRows` number of rows
            # note take() is an action
            # limit() is a transformation
            topDF = copyDF.limit( numRows )
            filePath = outputDir + "/" + str(i) + "/" + sampleName 
            self.writeCSV(df, logTag, filePath)

            
            # Truncate the `copy_df` to remove
            # the contents fetched for `temp_df`
            # original quant.sf files are sorted by name however
            # we must use order by, else the row names between
            # GTEx sample will not be the same
            # we can not simply sort or orderBy once. we have to 
            # do this on every iteration
            copyDF = copyDF.subtract(topDF).orderBy( "Name" )
            
            # Increment the split number
            i += 1
            
        if remainingRows > 0 :
            filePath = outputDir + "/" + str(i) + "/" + sampleName 
            self.writeCSV(copyDF, logTag, filePath)
                    
        self.logger.info( "END\n" )
        
    ###############################################################################
    # AEDWIP TODO rework, pass in list of files names
    # single for loop loads df and joins,
    def joinParts_deprecated(self, dfList):
        '''
        aedwip
        '''
        self.logger.info("BEGIN")
        retDF = dfList[0]

        for i in range( 1, len(dfList) ):
            df2 = dfList[i]
            retDF = retDF.join( df2.selectExpr("*"), on=["Name"] )
        
        # the columns should be sample names
        # it is important that they are in alphabetic order else
        # they will not match colData matrix used by DESeq2.
        # colData has meta data about each sample used as preditor variables 
        # by the glm
        #
        originalCols = retDF.columns
        self.logger.info("originalCols.sort() {}".format(originalCols.sort()))
        sortedCols = sorted( originalCols )

        retDF = retDF.select(sortedCols)
        self.logger.info("sorted retDF columns:{}".format(retDF.columns))

        self.logger.info("END\n")
        return retDF

    ###############################################################################
    def unionParts(self, dfList):
        '''
        aedwip
        '''
        self.logger.info("BEGIN")
        retDF = dfList[0]

        for i in range( 1, len(dfList) ):
            df2 = dfList[i]
            retDF = retDF.union(df2)

        self.logger.info("END\n")
        return retDF
