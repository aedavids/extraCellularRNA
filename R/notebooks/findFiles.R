#
# findFiles.R
# kras.ipsc DESeq utility functions
#
# top level functions
#   getKrasIpscDESeqDataSet <- function(rootDir, indexFilePath, design)

library("DESeq2")
library("tximport")

########################################################################
getKrasIpscDESeqDataSet <- function(rootDir,
                                    indexFilePath,
                                    findFilesFunc,
                                    getColDFFunc,
                                    design) {
  
  # arguments:
  #   rootDir: 
  #     root path to data file
  #
  #   indexFilePath:
  #     path to gene index file.
  #     example: /home/kimlab/genomes.annotations/gencode.32/gen.32.tx.2.gene.csv'
  #
  #   findFilesFunc
  #     example function: findDay7KrasIpscDataFiles()
  #
  #   getColDFFunc:
  #     example function: getKrasIPSCDay7ColDataDF()
  #
  # finds all the salmon quant files
  # checks to make sure salmon was run with same reference index file
  # uses tximport() to reads quant.sf files 
  # create column data frame (meta data describing each quant.sf file)
  #
  # return
  #   list$dds == DESeqDataSetFromTximport()
  #   list$fileList == data set files
  
  
  fileList <- findFilesFunc(rootDir)
  allSalmonQuantFiles <- unlist(fileList, use.names = FALSE)
  print( allSalmonQuantFiles )
  
  # make sure salmon was run with the sam reference index file
  isSameIndex(allSalmonQuantFiles)
  
  # AEDWIP parameterize this,
  transriptMapperDF <- transcriptMapper(indexFilePath, isHGNC=FALSE)
  
  # argument ignoreAfterBar whether to split the tx id on the '|' character to
  # facilitate matching with the tx id in transcriptMapperDF, ie tx2gene
  txi <- tximport(allSalmonQuantFiles, type="salmon", 
                  tx2gene=transriptMapperDF, ignoreAfterBar=TRUE)
  
  # get the meta data for each column of count data.
  colDF <- getColDFFunc(fileList)
  
  dds <- DESeqDataSetFromTximport(txi, colData=colDF, design)
  
  retList = List(dds, fileList)
  names(retList) <- c("dds", "fileList")
  return(retList) 
}

########################################################################
checkIndexIsSame <- function( cmdJSONList, refIdxFilePath) {
  # returns true if all the index files are the same
  
  idxIsSame = lapply(cmdJSONList, function(filePath){ filePath == refIdxFilePath })
  idxIsSame <- unlist( idxIsSame )
  return( all(idxIsSame) )
}

#
# function to find salmon quant.sf files
#
########################################################################
findConditionQuantFiles <- function (dirPath, condition="ctrl") {
  #dirPath /home/kimlab/kras.ipsc/data/bulk.data/day.5/
  # condition is typically either 'ctrl' or 'kras'
  # sample file
  # [1] "/home/kimlab/kras.ipsc/data/bulk.data/day.5/ctrl.1/gencode.salmon.out/quant.sf"
  
  f <- list.files( dirPath )
  conditionIdxList = grep(condition, f)
  conditionDirs = f[ conditionIdxList ]
  
  conditionSalmonFiles <- paste( dirPath, conditionDirs, 
                                     'gencode.salmon.out/quant.sf', 
                                     sep="/" )
  return( conditionSalmonFiles )
}

########################################################################
findDay5KrasIpscDataFiles <- function(rootDir) {
  #
  # arguments
  #   dirpath
  #     example: rootDir <- "/home/kimlab"
  #
  #
  # returns a list of file paths, name of list items will be dataSet name
  #
  
  dataFilesList = list()
  
  #### find all bulk day 5 control estimated transcript count files
  bulkDataDirPath = file.path( rootDir, 'kras.ipsc', 'data', 'bulk.data', 'day.5' )
  bulkDay5CntrlSalmonFiles <- findConditionQuantFiles(bulkDataDirPath, condition='ctrl')
  bulkDay5CntrlSalmonFiles
  if ( ! all(file.exists(bulkDay5CntrlSalmonFiles)) ) {
    stop(c("ERROR: file not found: ", bulkDay5CntrlSalmonFiles))
  }
  dataFilesList$bulkDay5CntrlSalmonFiles <- bulkDay5CntrlSalmonFiles
  
  #### find all day 5 exo control estimated transcript files
  exoDataDirPath = file.path( rootDir, 'kras.ipsc', 'data', 'exo.data', 'gen1c.day.5.exo.input.data' )
  exoDay5CntrlSalmonFiles <- findConditionQuantFiles(exoDataDirPath, condition='ctrl')
  exoDay5CntrlSalmonFiles
  if ( ! all(file.exists(exoDay5CntrlSalmonFiles)) ) {
    stop(c("ERROR: file not found: ", exoDay5CntrlSalmonFiles))
  }  
  dataFilesList$exoDay5CntrlSalmonFiles <- exoDay5CntrlSalmonFiles
  
  #### find all bulk day 5 kras estimated transcript count files
  bulkDay5KrasSalmonFiles <- findConditionQuantFiles(bulkDataDirPath, condition='kras')
  bulkDay5KrasSalmonFiles
  if ( ! all(file.exists(bulkDay5KrasSalmonFiles)) ) {
    stop(c("ERROR: file not found: ", bulkDay5KrasSalmonFiles))
  }
  dataFilesList$bulkDay5KrasSalmonFiles <- bulkDay5KrasSalmonFiles
  
  ### find all day 5 exo kras estimated transcript files
  exoDay5KrasSalmonFiles <- findConditionQuantFiles(exoDataDirPath, condition='kras')
  exoDay5KrasSalmonFiles
  if ( ! all(file.exists(exoDay5KrasSalmonFiles)) ) {
    stop(c("ERROR: file not found: ", exoDay5KrasSalmonFiles))
  }
  dataFilesList$exoDay5KrasSalmonFiles <- exoDay5KrasSalmonFiles
  
  return(dataFilesList)
}

########################################################################
findDay7KrasIpscDataFiles <- function(rootDir) {
  #
  # arguments
  #   dirpath
  #     example: rootDir <- "/home/kimlab"
  #
  #
  # returns a list of file paths, name of list items will be dataSet name
  #
  
  dataFilesList = list()
  
  #### find all bulk day 7 control estimated transcript count files
  bulkDataDirPath = file.path( rootDir, 'kras.ipsc', 'data', 'bulk.data', 'day.7' )
  bulkDay7CntrlSalmonFiles <- findConditionQuantFiles(bulkDataDirPath, condition='ctrl')
  bulkDay7CntrlSalmonFiles
  if ( ! all(file.exists(bulkDay7CntrlSalmonFiles)) ) {
    stop(c("ERROR: file not found: ", bulkDay7CntrlSalmonFiles))
  }
  dataFilesList$bulkDay7CntrlSalmonFiles <- bulkDay7CntrlSalmonFiles
  
  
  #### find all bulk day 7 kras estimated transcript count files
  bulkDay7KrasSalmonFiles <- findConditionQuantFiles(bulkDataDirPath, condition='kras')
  bulkDay7KrasSalmonFiles
  if ( ! all(file.exists(bulkDay7KrasSalmonFiles)) ) {
    stop(c("ERROR: file not found: ", bulkDay7KrasSalmonFiles))
  }
  dataFilesList$bulkDay7KrasSalmonFiles <- bulkDay7KrasSalmonFiles

  
  return(dataFilesList)
}

########################################################################
getKrasIPSCDay5ColDataDF <- function(fileList) {
  # create the colData data frame. This is the meta data describing each sample
  # The Rows of correspond to columns of countData and cols are features that
  # describe the replicates. The design model use the columns as features
  
  numBulkCntrl <- length( fileList$bulkDay5CntrlSalmonFiles )
  numExoCntrl  <- length( fileList$exoDay5CntrlSalmonFiles )
  numBulkKras  <- length( fileList$bulkDay5KrasSalmonFiles )
  numExoKras   <- length( fileList$exoDay5KrasSalmonFiles )
  
  # place holder for day 5 vs 7. we only have exo data for day 5 control
  day <- c( rep('5', numBulkCntrl + numExoCntrl + numBulkKras + numExoKras) ) 
  day <- factor( day )
  
  sampleType <- c( rep('bulk', numBulkCntrl), rep('exo', numExoCntrl), 
                   rep('bulk', numBulkKras), rep('exo', numExoKras))
  
  sampleType <- factor( sampleType )
  
  treatment <- c( rep('ctrl', numBulkCntrl + numExoCntrl),
                  rep('kras', numBulkKras + numExoKras)) 
  
  treatment <- factor(treatment)
  
  colDF <- data.frame(sampleType, treatment, day)
  
  # check levels. by default R choose alphabetic order
  # we need to make sure the control level is first or use the contrast argument
  levels(colDF$sampleType)
  levels(colDF$treatment)
  levels(colDF$day)
  
  return(colDF)  
}

########################################################################
getKrasIPSCDay7ColDataDF <- function(fileList) {
  # create the colData data frame. This is the meta data describing each sample
  # The Rows of correspond to columns of countData and cols are features that
  # describe the replicates. The design model use the columns as features
  
  numBulkCntrl <- length( fileList$bulkDay7CntrlSalmonFiles )
  numBulkKras  <- length( fileList$bulkDay7KrasSalmonFiles )
  
  
  # place holder for day 5 vs 7. we only have exo data for day 5 control
  day <- c( rep('7', numBulkCntrl + numBulkKras ) ) 
  day <- factor( day )
  
  sampleType <- c( rep('bulk', numBulkCntrl), rep('bulk', numBulkKras) )
  
  sampleType <- factor( sampleType )
  
  treatment <- c( rep('ctrl', numBulkCntrl ),
                  rep('kras', numBulkKras )) 
  
  treatment <- factor(treatment)
  
  colDF <- data.frame(sampleType, treatment, day)
  
  # check levels. by default R choose alphabetic order
  # we need to make sure the control level is first or use the contrast argument
  levels(colDF$sampleType)
  levels(colDF$treatment)
  levels(colDF$day)
  
  return(colDF)  
}

########################################################################
#
# utility functions to make sure all the
# salmon quant.sf files used the same index
#
getSalmonIndexFiles <- function(cmdInfoJsonFiles) {
  jResults <- lapply(cmdInfoJsonFiles, function(jsonFile){fromJSON(file=jsonFile)})
  idxList <- lapply(jResults, function(jr){ jr$index})
  
  flattenedList = unlist(idxList)
  return( flattenedList )
}

########################################################################
isSameIndex <- function(allSalmonQuantFiles) {
  # make sure salmon was run using the same index file 
  # calls stop if there is a mismatch
  
  # replace quant.sf with cmd_info.json
  cmdInfoJsonFiles <- lapply(allSalmonQuantFiles, function(quantFile) {
    paste( dirname( quantFile ), "cmd_info.json", sep="/" )
  } )
  
  # get a list of json objects
  cmdInfoJSONList <- getSalmonIndexFiles( cmdInfoJsonFiles )
  
  if ( ! checkIndexIsSame(cmdInfoJSONList, cmdInfoJSONList[1]) ) {
    stop( "ERROR: 'salmon quant' files used different index files!")
  } 
}

########################################################################
transcriptMapper <- function(indexFilePath, isHGNC=TRUE) {
  # returns a 2 column data frame that maps from TXTNAME to either HGNC or BIOTYPE
  # arguments:
  #.  indexFilePath:
  #     example:file.path( '/home/kimlab/genomes.annotations/gencode.32', 
  #                          'gen.32.tx.2.gene.csv')
  #   if isHGNC then maps transcript id to HGNC else BIOTYPE
  
  # tx2Gene.csv <- file.path( '/home/kimlab/genomes.annotations/gencode.32', 
  #                           'gen.32.tx.2.gene.csv')
  gen.32.tx.2.geneDF <- read.table(file=indexFilePath, header=FALSE, sep ="|")
  #print(head(gen.32.tx.2.geneDF))
  if (isHGNC) {
    #tx2geneDF <- gen.32.tx.2.geneDF[, c(1,5)] # 4 + 1 is HGNC with gene specific transcript id
    tx2geneDF <- gen.32.tx.2.geneDF[, c(1,6)] # 5 + 1 is HGNC required to stay at gene level
  } else  {
    tx2geneDF <- gen.32.tx.2.geneDF[, c(1,8)]  # 7 + 1 is bio type
  }
  
  names(tx2geneDF) <- c("TXNAME", "BIOTYPE")
  
  return( tx2geneDF )
}
