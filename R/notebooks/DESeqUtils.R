#
# DESeqUtils.R
# aedavids@ucsc.edu
#

library("DESeq2")
library("tximport")
source("findFiles.R") # assume file in in same directory as current Rmd

########################################################################
createComparableCounts <- function( dds, filterRowCount, filterSampleCount ) {
  # arguments
  #   dds: DESeq Data Set
  #
  # filterRowCount, filterSampleCount: integers
  #     filter rows arguments:
  #     remove rows if they do not have a count of filterRowCount or more in
  #     filterSampleCount or more samples
  #    . 
  # return normalized counts dds data frame
  #
  
  # warning("AEDWIP createComparableCounts() remove browser()", immediate. = TRUE)
  # browser()
  
  #warning("AEDWIP BEGIN createComparableCounts", immediate.=TRUE)
  
  # filter rows:
  # remove rows if they do not have a count of 5 or more in
  # 4 or more samples
  keep <- rowSums(counts(dds) >= filterRowCount) >= filterSampleCount
  table(keep)
  dds <- dds[keep,]  
  
  #warning("AEDWIP createComparableCounts completed filter", immediate.=TRUE)
  
  # warning("AEDWIP createComparableCounts remove browser()", immediate.=TRUE)
  # browser()
  
  # https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/estimateSizeFactors
  # This function estimates the size factors using the "median ratio method"
  # described by Equation 5 in Anders and Huber (2010). The estimated size
  # factors can be accessed using the accessor function sizeFactors
  dds <- estimateSizeFactors( dds )
  
  #warning("AEDWIP createComparableCounts completed estimateSizeFactors", immediate.=TRUE)
  
  
  # https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/counts
  #
  # The counts slot holds the count data as a matrix of non-negative integer
  # count values, one row for each observational unit (gene or the like), and
  # one column for each sample.
  #
  # normalized: logical indicating whether or not to divide the counts by the
  # size factors or normalization factors before returning (normalization
  # factors always preempt size factors)
  #
  # https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/normalizationFactors
  # Gene-specific normalization factors for each sample can be provided as a
  # matrix, which will preempt sizeFactors. In some experiments, counts for each
  # sample have varying dependence on covariates, e.g. on GC-content for
  # sequencing data run on different days, and in this case it makes sense to
  # provide gene-specific factors for each sample rather than a single size
  # factor.
  #
  #
  # https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/sizeFactors
  # The sizeFactors vector assigns to each column of the count matrix a value,
  # the size factor, such that count values in the columns can be brought to a
  # common scale by dividing by the corresponding size factor (as performed by
  # counts(dds, normalized=TRUE))
  
  ddsCounts <- counts(dds, normalized=TRUE)
  ddsCountsDF <- as.data.frame( ddsCounts )
  
  # warning("AEDWIP  createComparableCounts completed counts", immediate.=TRUE)
  # 
  # warning("AEDWIP END createComparableCounts", immediate.=TRUE)
  
  return( ddsCountsDF )
}

########################################################################
getDESeqDataSet <- function(rootDir,
                            indexFilePath,
                            findFilesFunc,
                            getColDFFunc,
                            design,
                            transcriptMapperDF,
                            ignoreAfterBar=TRUE,
                            txOut=FALSE) {
  
  # play around  with polymorphism by passing functions
  # TODO figure out how to deal with variable arguments.
  # example findA549_0.2MOI_24hrFiles(rootDir, index)
  
  #
  # ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta
  #
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
  #   design
  #     not used need to pass DESeq?()
  #
  #   transriptMapperDF:
  #     a data frame with 2 columns, first column is transcript Id.
  #     The second column is the id to map to. E.G. gene id, bio type, ...
  #     see
  #       findFiles.R genericTranscriptMapper()
  #       findFiles.R TE_TranscriptMapper()
  #       findFiles.R TE_bioTypeTranscriptMapper()
  #
  #   ignoreAfterBar=TRUE
  #       whether to split the tx id on the '|' character to facilitate matching 
  #       with the tx id in tx2gene
  #     
  #       if you are using gencode.v35.ucsc.rmsk.tx.to.gene.csv style mapping file
  #       pass FALSE, we want the entire name to be the transcript id
  #
  #   txOut=FALSE
  #       avoid gene-level summarization by setting txOut=TRUE
  #       if you are trying to normalize transcript and not map them to
  #       genes, or biotype, setting txOut=TRUE reduces required memory
  #       this allows use gencode.v35.ucsc.rmsk.tx.to.gene.csv

  # finds all the salmon quant files
  # checks to make sure salmon was run with same reference index file
  # uses tximport() to reads quant.sf files 
  # create column data frame (meta data describing each quant.sf file)
  #
  # return
  #   list$dds == DESeqDataSetFromTximport()
  #   list$fileList == data set files
  
  
  fileList <- findFilesFunc(rootDir)
  #print(sprintf("AEDWIP names(fileList) %s\n", names(fileList)))
  allSalmonQuantFiles <- unlist(fileList, use.names = FALSE)
  #print(sprintf("AEDWIP length(allSalmonQuantFiles) %s\n", length(allSalmonQuantFiles)))
  
  #print( allSalmonQuantFiles )
  
  # make sure salmon was run with the sam reference index file
  isSameIndex(allSalmonQuantFiles)
  
  #transriptMapperDF <- genericTranscriptMapper(indexFilePath, mapToId) 
  # transriptMapperDF <- transcriptMapper(indexFilePath, isHGNC=FALSE)
  
  
  
  # unable to run GUT data set using TE transcript mapping
  # work around is to set txOut=TRUE
  # AEDWIP TODO Roman WIP: instead of passing DF we can pass the csv file and will get more meta data
  # logical, whether the function should just output transcript-level (default FALSE)

  # tximport returns a list with matrices, "abundance", "counts", and "length", 
  # where the transcript level information is summarized to the gene-level. 
  # Typically, abundance is provided by the quantification tools as 
  # TPM (transcripts-per-million), while the counts are estimated counts 
  # (possibly fractional), and the "length" matrix contains the effective gene lengths
  #
  # argument ignoreAfterBar whether to split the tx id on the '|' character to
  # facilitate matching with the tx id in transcriptMapperDF, ie tx2gene
  #
  # avoid gene-level summarization by setting txOut=TRUE
  txi <- tximport(allSalmonQuantFiles, 
                  type="salmon", 
                  tx2gene=transcriptMapperDF, 
                  ignoreAfterBar=ignoreAfterBar,
                  txOut=txOut)


  
  # get the meta data for each column of count data.
  colDF <- getColDFFunc(fileList)
  
  dds <- DESeqDataSetFromTximport(txi, colData=colDF, design)
  
  retList = List(dds, fileList)
  names(retList) <- c("dds", "fileList")
  
  return(retList) 
}

########################################################################
loadAndNormlize <- function(parameters) {
  # 
  # return a list ddsCountsDF", "dds", "fileList"
  # retList$ddsCountsDF <- counts data frame
  # retList$dds <- DESeqData object
  # retList$fileList <- list list of salmon quant files, inner list is partitioned by condition
  #
  
  
  
  # example mapToId <- SALMON_QUANT_SF_NAMES_LIST()$ENSG
  # mapToId <- SALMON_QUANT_SF_NAMES_LIST()$ENSG
  # mapToIdName <- names(SALMON_QUANT_SF_NAMES_LIST()[mapToId])
  
  # mapToIdName <- "tx"
  # transcriptMapperDF <- TE_TranscriptMapper(tx2MappingFilePath)
  # 
  # ignoreAfterBar <- False
  # warning we are not running,using the model. 
  
  # debug
  print("loadAndNormlize input parameters")
  print(parameters)
  
  #
  # to all the config and directory creation first
  # you do not want to run for a couple of hours only to have a config
  # bug prevent you from saving results
  #
  
  tx2MappingFile <- parameters$tx2MappingFile
  # outputDir <- paste0('normalizedCounts.', tx2MappingFile)
  # #outputRoot <- file.path("/home/aedavids/extraCellularRNA/data/R/output", outputDir)
  # outputRoot <- parameters$outputRoot
  # outputRoot <- file.path(outputRoot, outputDir)
  # dir.create(outputRoot, recursive=TRUE, showWarnings = FALSE)
  
  #tx2MappingFilePath <- file.path(parameters$tx2MappingDir, parameters$tx2MappingFile)
  # mapToIdName <- parameters$mapToIdName
  ignoreAfterBar <- parameters$ignoreAfterBar
  
  tx2MappingDir <- parameters$tx2MappingDir
  # print(names(parameters))
  # warning(sprintf("tx2MappingDir:%s tx2MappingFile:%s\n", tx2MappingDir, tx2MappingFile),
  #         immediate.=TRUE)
  # rsucks <- is.function(parameters$transcriptMapperDFFunc)
  # warning(sprintf("AEDWIP is_function: %d\n", rsucks ), immediate. = TRUE)
  # warning(parameters$transcriptMapperDFFunc, immediate.=TRUE)
  # warning("aEDWIP r sucks remove browser\n")
  # browser()
  transcriptMapperDF <- parameters$transcriptMapperDFFunc(tx2MappingDir, tx2MappingFile)
  
  #csvOutFileName <- "pancreas.plasma.ev.long.RNA.normalized.deseq.biotype.counts.csv"
  # csvOutFileName <- paste0("pancreas.plasma.ev.long.RNA.normalized.deseq.", 
  #                          mapToIdName, 
  #                          ".counts.csv")  
  
  #rootDir <- parameters$rootDir
  
  #warning("loadAndNormlize() before getDESeqDataSet\n", immediate.=TRUE) 
  # warning("AEDWIP TODO loadAndNormlize() parameterize findFilesFunc and getColDFFunc\n", 
  #         immediate.=TRUE)
  #warning("AEDWIP TODO loadAndNormlize() rootDir is a global!!!!\n", 
  
  # warning("aedwip remove browser load()")
  # browser()
  ddsList <- getDESeqDataSet(rootDir=parameters$rootDir,
                             tx2MappingFilePath,
                             findFilesFunc= parameters$findFilesFunc, 
                             getColDFFunc=parameters$getColDFFunc,
                             design=~ diseaseState,
                             transcriptMapperDF,
                             ignoreAfterBar=ignoreAfterBar,
                             txOut=parameters$txOut)
  
  # warning("AEDWIP remove browser()", immediate. = TRUE)
  # browser()
  
  dds <- ddsList$dds
  fileList <- ddsList$fileList
  
  #warning("loadAndNormlize() after getDESeqDataSet\n", immediate.=TRUE)  
  
  #
  # pre filtering:
  # ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
  # pre-filter low count genes is not necessary.
  # the reasons which make pre-filtering useful: 
  #   1) reduce the memory of dds data object
  #   2) increase the speed of the transformation and testing functions within DESeq2
  #
  # https://support.bioconductor.org/p/65256/#65260
  # "you typically don't need to pre-filter because independent filtering occurs
  # within results() to save you from multiple test correction on genes with no
  # power (see ?results and the vignette section about independent filtering, or
  # the paper). The main reason to pre-filter would be to increase speed.
  # Designs with many samples and many interaction terms can slow down on genes
  # which have very few reads.
  #
  # reson results() filters is having fewer genes improves the adjusted p-values.

  
  #aediwp 5,4 was when we had 6 replicants in 2 conditions
  print("AEDWIP TODO check filtering see statquest https://youtu.be/Gi0JdrxRq5s edgeR and DESeq2, part2: Independent Filtering (removing genes with low read counts)" )
  
  # minNumberOfSamples(aedwip)
  # df <- createComparableCounts( dds, filterRowCount=5, filterSampleCount=4 )
  
  filterSampleCount = minNumberOfSamples(fileList, percent=0.25)
  print(sprintf("filterSampleCount:%f\n", filterSampleCount))
  aedwip <- 5 #with 5 we got over 17k genes using the only 4 quant.sf files
  #warning("normalizeRunner before createComparableCounts", immediate.=TRUE)  
  ddsCountsDF <- createComparableCounts( dds, filterRowCount=aedwip, filterSampleCount )
  ddsCountsDF <- createComparableCounts( dds, filterRowCount=aedwip, filterSampleCount )
  #warning("normalizeRunner after createComparableCounts", immediate.=TRUE)  
  
  
  warning("AEDWIP hack to aid debuging, ")
  ret = list(ddsCountsDF,  dds, fileList) 
  names(ret) <- c("ddsCountsDF", "dds", "fileList")
  
  return(ret)
}

########################################################################
saveNormalizedDF <- function(dsCountsDF, parameters) {
  # broke this out to make debugging easier
  
  tx2MappingFile <- parameters$tx2MappingFile
  outputDir <- paste0('normalizedCounts.', tx2MappingFile)
  #outputRoot <- file.path("/home/aedavids/extraCellularRNA/data/R/output", outputDir)
  outputRoot <- parameters$outputRoot
  outputRoot <- file.path(outputRoot, outputDir)
  dir.create(outputRoot, recursive=TRUE, showWarnings = FALSE)
  
  mapToIdName <- parameters$mapToIdName
  
  stop("AEDWIP TODO paramterize out file name it is hard code")
  csvOutFileName <- paste0("pancreas.plasma.ev.long.RNA.normalized.deseq.", 
                           mapToIdName, 
                           ".counts.csv")   
  
  csvFilePath <- file.path( outputRoot, csvOutFileName)
  # use cat instead of print if you want to add new lines
  # https://stackoverflow.com/a/9317914/4586180
  # print(sprintf("save to:\n%s\n",csvFilePath))
  cat(sprintf("save to:\n%s\n",csvFilePath))
  
  # write_csv does not write row names
  write.csv(dsCountsDF, csvFilePath)
}
