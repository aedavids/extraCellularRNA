#!/usr/bin/env Rscript
# Distributed DESeq2 Test as script
# test deseq works from script. Uses a small hacked test data set
# need to execute a hack work around to run DESeq. 
# ref: extraCellularRNA/R/notebooks/reverseEngineerDESeq/distrubutedDESeq2Test.Rmd

# capture output file
sink('testDESeqScript.out')

library("argparse")
library("DESeq2")
library("BiocParallel")
#
# parse cli
# https://cran.r-project.org/web/packages/argparse/vignettes/argparse.html
#
parser <- ArgumentParser(description='AEDWIP test docker')

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


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args() # c('--countMatrix', 'aedwipFile'))
cat("\nargs\n")
print(args)

countMatrixFile <- args$countMatrix
cat("\ncountMatrixFile \n")
print(countMatrixFile)

colDataFile <- args$colData
cat("\ncolDataFile\n")
print(colDataFile)

design <- args$design
cat("\ndesign \n")
print(design)

referenceLevel <- args$referenceLevel
cat("\nreferenceLevel \n")
print(referenceLevel)

outFile <- args$outFile
cat("\noutFile \n")
print(outFile)

numCores <-args$numCores
cat("\nnumCores\n")
print(numCores)

#
# find the variable of interest in the design formula
# it is the last variable in the formula. 
# strange syntax to parse on whitespace
#
tokens <- strsplit(design, " +")[[1]]
variableOfInterest <-  tokens[length(tokens)]
cat("\n \n")
print(variableOfInterest)

#
# load the count Matrix data file
#
if( file.access(countMatrixFile) == -1) {
     stop(sprintf("Specified  count matrix file ( %s ) does not exist", countMatrixFile))
} else {
  countMatrixDF <- read.table(countMatrixFile, sep="\t", header=TRUE)
  cat("\n \n")
  print(countMatrixDF)
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
  
  cat("\ncolData \n")
  print( colDataDF )
}


# convert the data frame to matrix
countMatrix <- data.matrix( countMatrixDF )
cat("\ncount matrix\n")
print(countMatrix)



#
# run DESeq
#

#
# config DESeq to run in parallel
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#speed-up-and-parallelization-thoughts
#
# add parameter parallel=TRUE to DESeq(), results(), and lfcShrink()
register(MulticoreParam(numCores))

cat( "\nStep 1, run through the vignette quick start example\n" )


# create a DESeqDataSetFromMatrix object
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = colDataDF,
                              design = as.formula(design) )



cat( "\nraw counts\n")
counts(dds)

cat( "\ndesign\n")
design( dds )


# our test data is funky
cat("\nhack work around for or funky test data\n")

dds <- estimateSizeFactors(dds)

cat("\n sizeFactors\n")
sizeFactors( dds )

cat("\n now that we have size factors we can get normalized counts\n")

print( counts(dds, normalized=TRUE) ) # get()

cat("\nwork around step 2) Estimate the dispersions for a DESeqDataSet work around\n")
dds <- estimateDispersionsGeneEst(dds)
print( mcols(dds) )

dispersions(dds) <- mcols(dds)$dispGeneEst # set
cat("\n gene dispersions\n")
dispersions(dds)


cat("\nwork around step 3) Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest\n" )
dds <- nbinomWaldTest(dds)
cat("\nmcols\n")
print( mcols(dds) )

# get a list of result coefficients
cat("\nresultsName\n")
print( resultsNames( dds ) )



cat("\nsave our expected results\n")
# results takes a name argument. by default we sampleType treatment vs control
# see the documentation there are lots of arguments
expectedDESeqResults <- results(dds, parallel=TRUE)
expectedDESeqResults

write.table( as.data.frame(expectedDESeqResults),
           file = outFile,
           sep="\t" )

# turn output capture off
sink()

# exit session do not save
q(save="no")

