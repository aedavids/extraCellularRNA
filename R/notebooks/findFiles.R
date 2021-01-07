#
# findFiles.R
# typically used by DESeq utility functions
#
# # top level functions
# findDay5KrasIpscDataFiles <- function(rootDir)
# findDay7KrasIpscDataFiles <- function(rootDir)
# findPancreasPlasmaEV_longRNAFiles <- function(rootDir)
# genericTranscriptMapper (indexFilePath, SALMON_QUANT_SF_NAMES_LIST_item)
# transcriptMapper(indexFilePath, isHGNC=TRUE) 
#
# global:
# https://github.com/rreggiar/welcome-to-the-kim-lab/wiki/Salmon-output-data-dictionary#understanding-the-name-column
# SALMON_QUANT_SF_NAMES_LIST

library('dplyr')
library('readr') # tximport() will use readr if present
library('rjson')
library('tidyr')
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
getPancreasPlasmaEV_longRNAColDataDF <- function(fileList) {
  # create the colData data frame. This is the meta data describing each sample
  # The Rows of correspond to columns of countData and cols are features that
  # describe the replicates. The design model use the columns as features
  
  numHealthy <- length( fileList$healthySalmonFiles )
  numPDAC  <- length( fileList$PDACSalmonFiles )
  
  diseaseState <- c( rep('healthy', numHealthy ),
                     rep('PDAC', numPDAC )) 
  
  diseaseState <- factor(diseaseState, levels=c("healthy", "PDAC"))
  
  colDF <- data.frame(diseaseState)
  
  return(colDF)  
}

########################################################################
findPancreasPlasmaEV_longRNAFiles <- function(rootDir) {
  #
  # finds salmon quant.sf files. Checks to make sure file exists
  #
  # arguments
  #   dirpath
  #     example: rootDir <- "/home/kimlab/pancreas.plasma.ev.long.RNA"
  #
  #
  # returns a list of file paths, name of list items will be dataSet name
  #   with names healthySalmonFiles, PDACSalmonFiles
  #
  
  # warning("AEDWIP findPancreasPlasmaEV_longRNAFiles remove browser")
  # browser()
  
  dataFilesList = list()
  
  dataFilesList$healthySalmonFiles <- findPancreasPlasmaEVForCondition(rootDir, "healthy")
  dataFilesList$PDACSalmonFiles <- findPancreasPlasmaEVForCondition(rootDir, "PDAC")
  
  # warning("AEDWIP DEBUG findPancreasPlasmaEV_longRNAFiles() return subset of files to make debug faster\n",
  #         immediate.=TRUE)
  # dataFilesList$healthySalmonFiles <- findPancreasPlasmaEVForCondition(rootDir, "healthy")[1:2]
  # dataFilesList$PDACSalmonFiles <- findPancreasPlasmaEVForCondition(rootDir, "PDAC")[1:2]
  
  print(sprintf("length(dataFilesList$healthySalmonFiles) = %d", length(dataFilesList$healthySalmonFiles)))
  print(sprintf("length(dataFilesList$PDACSalmonFiles) =  %d", length(dataFilesList$PDACSalmonFiles)))
  
  # warning("AEDWIP DEBUG findPancreasPlasmaEV_longRNAFiles()\n")
  # print(sprintf("healthySalmonFiles\n%s", 
  #               dataFilesList$healthySalmonFiles))
  # print(sprintf("PDACSalmonFiles\n%s", 
  #               dataFilesList$PDACSalmonFiles))
  
  return( dataFilesList)
}

########################################################################
findPancreasPlasmaEVForCondition <- function(rootDir, condition) {
  #
  # finds salmon quant.sf files. Checks to make sure file exists
  #
  # arguments
  #   rootDir: example  "/home/kimlab/pancreas.plasma.ev.long.RNA"
  #
  #   condition: "PDAC" or "healthy"
  #
  # returns
  #   a list of full salmon file path
  dataDirPath <- file.path( rootDir, "data", condition)
  sampleIds <- list.files( dataDirPath )
  
  salmonFiles <- file.path( dataDirPath, sampleIds, "salmon.out", "quant.sf" )
  if ( ! all(file.exists(salmonFiles)) ) {
    stop(c("ERROR: file not found: ", salmonFiles))
  }
  
  return(salmonFiles)
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
minNumberOfSamples = function(fileList, percent=0.25) {
  # returns filterSampleCount value for use with createComparableCounts()
  # use so that we only consider transcript that are in at % of our samples
  #
  # arguments:
  #   percent, a float > 0 and < 1
  #
  # see statquest https://youtu.be/Gi0JdrxRq5s edgeR and DESeq2, 
  # part2: Independent Filtering (removing genes with low read counts)
  

  fileCount <- 0
  for (lname in names(fileList)) {
    fileCount <- fileCount + length(fileList[[lname]])
  }
  fileCount <- as.integer(round( fileCount * percent ))
  
  return( fileCount )
}

########################################################################
# AEDWIP
# https://github.com/rreggiar/welcome-to-the-kim-lab/wiki/Salmon-output-data-dictionary#understanding-the-name-column
SALMON_QUANT_SF_NAMES_LIST <- function() {
  # https://github.com/rreggiar/welcome-to-the-kim-lab/wiki/Salmon-output-data-dictionary#understanding-the-name-column
  
  ret <- list()
  ret$ENST       <- 1
  ret$ENSG       <- 2
  ret$OTTHUMG    <- 3
  ret$OTTHUMT    <- 4
  ret$HGNCT      <- 5
  ret$HGNC       <- 6
  ret$REF_LENGTH <- 7   
  ret$BIO_TYPE   <- 8  
  
  return( ret )
}

########################################################################
getBiotypeMapper <- function(tx2MappingDir, tx2MappingFile) {
  # function passes in parameter list
  # enable lazy evaluation
  tx2MappingFilePath <- file.path(tx2MappingDir, tx2MappingFile)
  
  mapToId <- unlist( SALMON_QUANT_SF_NAMES_LIST()["BIO_TYPE"] )
  ret <- genericTranscriptMapper(tx2MappingFilePath, mapToId)
  
  return (ret)
}

########################################################################
getENSTMapper <- function(tx2MappingDir, tx2MappingFile) {
  # function passes in parameter list
  # enable lazy evaluation
  tx2MappingFilePath <- file.path(tx2MappingDir, tx2MappingFile)
  
  mapToId <- unlist( SALMON_QUANT_SF_NAMES_LIST()["ENST"] )
  ret <- genericTranscriptMapper(tx2MappingFilePath, mapToId)

  return (ret)
}

########################################################################
getHGNCMapper <- function(tx2MappingDir, tx2MappingFile) {
  # function passes in parameter list
  # enable lazy evaluation
  tx2MappingFilePath <- file.path(tx2MappingDir, tx2MappingFile)
  
  mapToId <- unlist( SALMON_QUANT_SF_NAMES_LIST()["HGNC"] )
  ret <- genericTranscriptMapper(tx2MappingFilePath, mapToId)
  
  return (ret)
}


########################################################################
genericTranscriptMapper <- function(indexFilePath, mapToId) {
  # returns a 2 column data frame that maps from TXTNAME to either HGNC or BIOTYPE
  # arguments:
  #.  indexFilePath:
  #     example:file.path( '/home/kimlab/genomes.annotations/gencode.32', 
  #                          'gen.32.tx.2.gene.csv')
  #   mapToId
  #     an item from the SALMON_QUANT_SF_NAMES_LIST
  #     example SALMON_QUANT_SF_NAMES_LIST$BIO_TYPE
  #     this id will be the second column in return data frame
  
  # tx2Gene.csv <- file.path( '/home/kimlab/genomes.annotations/gencode.32', 
  #                           'gen.32.tx.2.gene.csv')
  
  # sample data from
  # kimlab/genomes.annotations/gencode.35/gencode.v35.tx.to.gene.csv
  # ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1
  #   |DDX11L1-202|DDX11L1|1657|processed_transcript|,DDX11L1
  # ENST00000450305.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000002844.2
  #   |DDX11L1-201|DDX11L1|632|transcribed_unprocessed_pseudogene|,DDX11L1
  
  # readr read_table is faster
  gen.32.tx.2.geneDF <- read.table(file=indexFilePath, header=FALSE, sep ="|")
  tx2geneDF <- gen.32.tx.2.geneDF[, c(1,mapToId)]

  #print(head(gen.32.tx.2.geneDF))
  
  n <- names(SALMON_QUANT_SF_NAMES_LIST()[mapToId])
  names(tx2geneDF) <- c("TXNAME", n)
  
  return( tx2geneDF )
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
  
  # readr read_table is faster
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

