#!/usr/bin/env Rscript
# Distributed DESeq2 Test as script
# test deseq works from script. Uses a small hacked test data set
# need to execute a hack work around to run DESeq. 
# ref: extraCellularRNA/R/notebooks/reverseEngineerDESeq/distrubutedDESeq2Test.Rmd

# capture output file
sink('DESeqScript.out')

library("argparse")
library("DESeq2")
library("BiocParallel")
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
} else {
  countMatrixDF <- read.table(countMatrixFile, sep="\t", header=TRUE)
  cat("\n head(countMatrixDF)\n")
  print( head(countMatrixDF) )
}

# load the colData file
if( file.access(colDataFile) == -1) {
  stop(sprintf("Specified  colData ( %s ) does not exist", colDataFile))
} else {
  colDataDF <- read.table(colDataFile, 
                          sep="\t", 
                          header=TRUE,
                          stringsAsFactors = TRUE
                          )
  
  # by default R orders factors alphabetically. This will mess up our results
  cat( "\nlevels before relevel\n")
  print( levels( colDataDF[variableOfInterest][,1] ) )
  
  colDataDF[variableOfInterest][,1] <- relevel(colDataDF[variableOfInterest][,1], referenceLevel)
  
  cat("\nnew levels\n")
  print( levels(colDataDF[variableOfInterest][,1]) )
  
  cat("\nhead(colData) \n")
  print( head(colDataDF) )
  
  # https://github.com/mikelove/tximport/blob/master/R/tximport.R
  # line 361
}


# convert the data frame to matrix
# [,-1] means drop the first column
countMatrix <- data.matrix( countMatrixDF[,-1] )
cat("\n head(countMatrix) \n")
print( head(countMatrix) )

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



cat( "\nraw head(counts)\n")
head( counts(dds) )

dds <- estimateSizeFactors(dds)

cat("\n Step 2, get estimate size Factors\n")
write.table( as.data.frame( sizeFactors( dds ) ),
             file = estimateSizeFactorsOutfile,
             sep="\t" )

cat("\n Step 3,  head(normalized counts)\n")

print( head(counts(dds, normalized=TRUE)) ) # get()

# cat("\nwork around step 2) Estimate the dispersions for a DESeqDataSet work around\n")
# dds <- estimateDispersionsGeneEst(dds)
# print( mcols(dds) )
# 
# dispersions(dds) <- mcols(dds)$dispGeneEst # set
# cat("\n gene dispersions\n")
# dispersions(dds)

cat("\n Step 4, estimate dispersions")
dispersions(dds) <- estimateDispersions( dds )
cat("\n head(mcols(dds))\n")
print( head(mcols(dds)) )


cat("\n Step 5,  Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest\n" )
dds <- nbinomWaldTest(dds)
cat("\nmcols\n")
print( mcols(dds) )

# get a list of result coefficients
cat("\nresultsName\n")
print( resultsNames( dds ) )



cat("\nsave our expected results\n")
# results takes a name argument. by default we sampleType treatment vs control
# see the documentation there are lots of arguments
DESeqResults <- results(dds, parallel=TRUE)
DESeqResults

write.table( as.data.frame(DESeqResults),
           file = outFile,
           sep="\t" )

# turn output capture off
sink()

# exit session do not save
q(save="no")

