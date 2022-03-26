#!/usr/bin/env Rscript
# Distributed DESeq2 Test as script
# test deseq works from script. Uses a small hacked test data set
# need to execute a hack work around to run DESeq. 
# ref: extraCellularRNA/R/notebooks/reverseEngineerDESeq/distrubutedDESeq2Test.Rmd

# capture output file
sink('DESeqScript.out')

startTime <- Sys.time()

library("argparse")
library("DESeq2")
library("BiocParallel")

saveResults <- function(outFile, dds, DESeqResults ) {
  # write self describing meta data to header section of output file
  descriptionList <- DESeqResults@elementMetadata@listData[["description"]]
  cat( sprintf("%s \n", descriptionList[[1]]), file=outFile)
  for (i in 2:length(descriptionList)) {
      txt <- sprintf( "%s \n", descriptionList[[i]] ) 
      cat( txt, file=outFile, append=TRUE)
  }

  txt <- sprintf( "design: %s \n", format(design(dds)))
  cat( txt, file = outFile, append = TRUE)

  # create a column named 'name' with the gene name
  # else the output from write table will header line will not have
  # the same number of columns as the data rows

  cat("debug rownames(DESeqResults) \n")
  # cat( rownames(DESeqResults) )
  # cat(" \n")

  DESeqResults<- cbind( rownames(DESeqResults), DESeqResults )
  cat("debug  rownames after cbind( rownames(DESeqResults), DESeqResults )\n")
  # cat( rownames(DESeqResults) )

  colnames(DESeqResults)[1] <- "name"
  cat("debug  rownames setting first col to 'name'\n")
  #cat( rownames(DESeqResults) )

  cat("debug before write.table head(DESeqResults)\n")
  print( head(DESeqResults) )

  cat("debug callign write.table\n")
  write.table( DESeqResults,
               file = outFile,
               sep=delimator,
               row.names=FALSE,
               append=TRUE)
}

#
# parse cli
# https://cran.r-project.org/web/packages/argparse/vignettes/argparse.html
#
parser <- ArgumentParser(description='Runs DESeq2 differnetial experssion')

parser$add_argument("-c", "--countMatrix", 
                    required=TRUE,
                    help="count matrix in tsv format, columns headers are sample names, first row is gene name")

parser$add_argument("-m", "--colData", 
                    required=TRUE,
                    help="colData in tsv format, meta data about each sample")

parser$add_argument("-d", "--design",
                    required=TRUE,
                    help="design formula. All strings are assumed to be factors. Put the variable of interestes last. Example, '~ age + treatment'" )

parser$add_argument("-r", "--referenceLevel",
                    required=TRUE,
                    help="to compare C vs B, make B the reference level" )

parser$add_argument("-o", "--outFile",
                    required=TRUE,
                    help="file to write results to" )

parser$add_argument("-n", "--numCores",
                    type="integer",
                    default=1,
                    help="set > 1 to enable child process parallization" )

parser$add_argument("-e", "--estimateSizeFactorsOutfile", 
                    required=TRUE,
                    help="file to write estimated size factors used to scale counts" )

parser$add_argument("-1", "--oneVsAll",
                    required=FALSE,
                    action='store_true',
                    default=FALSE,
                    help="test reference Level vs. not reference Level")

parser$add_argument("-i", "--isCSV",
                    required=FALSE,
                    action='store_true',
                    default=FALSE,
                    help="default assume tsv files, if flag argument is present assume tcv files")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args() # c('--countMatrix', 'aedwipFile'))
cat("\nargs\n")
print(args)

countMatrixFile <- args$countMatrix
colDataFile <- args$colData
design <- args$design
referenceLevel <- args$referenceLevel
outFile <- args$outFile
numCores <-args$numCores
estimateSizeFactorsOutfile <- args$estimateSizeFactorsOutfile
oneVsAll <- args$oneVsAll

if (args$isCSV) {
  delimator <- ","
} else {
  delimator <- "\t"
}

#
# find the variable of interest in the design formula
# it is the last variable in the formula. 
# strange syntax to parse on whitespace
#
tokens <- strsplit(design, " +")[[1]]
variableOfInterest <-  tokens[length(tokens)]
cat( sprintf("\n variableOfInterest: %s \n", variableOfInterest) )

#
# load the count Matrix data file
#
if( file.access(countMatrixFile) == -1) {
     stop(sprintf("Specified  count matrix file ( %s ) does not exist", countMatrixFile))
} 

countMatrixDF <- read.table(countMatrixFile, sep=delimator, header=TRUE)
cat("\n head(countMatrixDF)[,1:3]\n")
print( head(countMatrixDF)[,1:3] )


# load the colData file
if( file.access(colDataFile) == -1) {
  stop(sprintf("Specified  colData ( %s ) does not exist", colDataFile))
} 

colDataDF <- read.table(colDataFile, 
                        sep=delimator, 
                        header=TRUE,
                        stringsAsFactors = TRUE
                        )

cat("\n debug original colData meta data")
#print( str(colData) )

# by default R orders factors alphabetically. This will mess up our results
cat( "\ndebug levels before relevel\n")
#print( levels( colDataDF[variableOfInterest][,1] ) )

colDataDF[variableOfInterest][,1] <- relevel(colDataDF[variableOfInterest][,1], referenceLevel)

cat("\n debug new levels\n")
# print( levels(colDataDF[variableOfInterest][,1]) )

if (oneVsAll) {
  notLevel = paste0("not_", referenceLevel)
  cat("\ndebug not level")
  print(notLevel)
  #levels( colDataDF[variableOfInterest] )[ levels(colDataDF[variableOfInterest] ) != referenceLevel ] <- notLevel
  # https://stackoverflow.com/questions/19730806/access-data-frame-column-using-variable
  levels( colDataDF[, variableOfInterest] )[ levels(colDataDF[, variableOfInterest] ) != referenceLevel ] = notLevel

  colDataDF[variableOfInterest][,1] <- relevel(colDataDF[variableOfInterest][,1], notLevel)


  cat("\none vs. all levels\n")
  print( levels(colDataDF[variableOfInterest][,1] ) )
  print("str(colDataDF)")
  print(str(colDataDF))
  print("colDataDF")
  print(colDataDF)
  print("####### end 1 vs. all debug ##########")
}

cat("\nhead(colData) \n")
print( head(colDataDF) )

# https://github.com/mikelove/tximport/blob/master/R/tximport.R
# line 361


# convert the data frame to matrix
# [,-1] means drop the first column
countMatrix <- data.matrix( countMatrixDF[,-1] )

# add row names to matrix so that gene names will be printed in final output
# row names are not considered data in the actual matrix
geneNameVector <- countMatrixDF$geneId
rownames( countMatrix ) <- geneNameVector

cat("\n head(countMatrix)[,1:3] \n")
print( head(countMatrix)[,1:3] )

#
# run DESeq
#

#
# config DESeq to run in parallel
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#speed-up-and-parallelization-thoughts
#
# add parameter parallel=TRUE to DESeq(), results(), and lfcShrink()
register(MulticoreParam(numCores))

cat( "\nStep 1, Create a DESeqDataSetFromMatrix\n" )


# create a DESeqDataSetFromMatrix object
# use round to convert numReads from double to integer
dds <- DESeqDataSetFromMatrix(countData = round(countMatrix),
                              colData = colDataDF,
                              design = as.formula(design) )



cat( "\nraw head(counts)[,1:3]\n")
head( counts(dds) )[,1:3]

dds <- estimateSizeFactors(dds)

cat("\n Step 2, get estimate size Factors\n")
write.table( as.data.frame( sizeFactors( dds ) ),
             file = estimateSizeFactorsOutfile,
             sep=delimator,
             row.names=FALSE )

cat("\n Step 3,  head(normalized counts)[,1:3]\n")

print( head(counts(dds, normalized=TRUE))[,1:3] ) # get()

# cat("\nwork around step 2) Estimate the dispersions for a DESeqDataSet work around\n")
# dds <- estimateDispersionsGeneEst(dds)
# print( mcols(dds) )
# 
# dispersions(dds) <- mcols(dds)$dispGeneEst # set
# cat("\n gene dispersions\n")
# dispersions(dds)

cat("\n Step 4, estimate dispersions")
dds <- estimateDispersions( dds )
cat("\n head(mcols(dds))[,1:3]\n")
print( head(mcols(dds)) [,1:3])


cat("\n Step 5,  Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest\n" )
dds <- nbinomWaldTest(dds)
cat("\n head(mcols))\n")
print( head(mcols(dds)) )

# get a list of result coefficients
cat("\nresultsName\n")
print( resultsNames( dds ) )



cat("\nsave our expected results\n")
# results takes a name argument. by default we sampleType treatment vs control
# see the documentation there are lots of arguments
DESeqResults <- results(dds, parallel=TRUE)

# sort by adjusted p-value
DESeqResults <- DESeqResults[order(DESeqResults$padj),]

print( head(DESeqResults) )



# # write self describing meta data to header section of output file
# descriptionList <- DESeqResults@elementMetadata@listData[["description"]]
# cat( sprintf("%s \n", descriptionList[[1]]), file=outFile)
# for (i in 2:length(descriptionList)) {
#   txt <- sprintf( "%s \n", descriptionList[[i]] ) 
#   cat( txt, file=outFile, append=TRUE)
# }

# txt <- sprintf( "design: %s \n", format(design(dds)))
# cat( txt, file = outFile, append = TRUE)

# # create a column named 'name' with the gene name
# # else the output from write table will header line will not have
# # the same number of columns as the data rows

# cat("debug rownames(DESeqResults) \n")
# # cat( rownames(DESeqResults) )
# # cat(" \n")

# DESeqResults<- cbind( rownames(DESeqResults), DESeqResults )
# cat("debug  rownames after cbind( rownames(DESeqResults), DESeqResults )\n")
# # cat( rownames(DESeqResults) )

# colnames(DESeqResults)[1] <- "name"
# cat("debug  rownames setting first col to 'name'\n")
# #cat( rownames(DESeqResults) )

# cat("debug before write.table head(DESeqResults)\n")
# print( head(DESeqResults) )

# cat("debug callign write.table\n")
# write.table( DESeqResults,
#            file = outFile,
#            sep=delimator,
#            row.names=FALSE,
#            append=TRUE)

saveResults( outFile, dds, DESeqResults )

Cat("DESeqScript.R completed sucessfully\n")

# turn output capture off
sink()

endTime <- Sys.time()
cat("\nrun time")
print( endTime - startTime )

# exit session do not save
q(save="no")

