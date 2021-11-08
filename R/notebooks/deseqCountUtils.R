#
# deseqCountUtils.R
#
# aedavids@ucsc.edu
# 6/1/21
#
source( "findFiles.R" ) 
source( "te.R" ) 

debugParametersList <- list(
  'debug'=list(
    'rootDir' = "/home/aedavids/extraCellularRNA/data/debugTximport/pancreas.plasma.ev.long.RNA",
    'findFilesFunc' = findPancreasPlasmaEV_longRNAFiles, 
    'getColDFFunc' = getPancreasPlasmaEV_longRNAColDataDF,    
    'tx2MappingFile' = 'debugTx2Gene.csv',
    'tx2MappingDir' = file.path('/home/aedavids/extraCellularRNA/data/debugTximport',
                                'genomes.annotations/gencode.v35'),  
    # 'outputRoot' = "/home/aedavids/extraCellularRNA/data/R/output/debug",
    'outputRoot' = "/home/aedavids/extraCellularRNA/data/R/debugTximport/output",
    'mapToIdName'='gene',
    'normalizedFileName'='pancreas.plasma.ev.long.RNA.normalized.deseq',
    'transcriptMapperDFFunc'=getTEGeneMapper,    
    
    # ignoreAfterBar logical, whether to split the tx id on the '|' character to 
    # facilitate matching with the tx id in tx2gene (default FALSE)
    #'ignoreAfterBar'=FALSE #bug we do not get gene names in output
    'ignoreAfterBar'=FALSE,
    
    #  txOut     avoid gene-level summarization by setting txOut=TRUE
    #       if you are trying to normalize transcript and not map them to
    #       genes, or biotype, setting txOut=TRUE reduces required memory
    #       this allows use gencode.v35.ucsc.rmsk.tx.to.gene.csv    
    'txOut'=FALSE,
    'designStr' = "~ diseaseState"
  ),
  'TE.bioType'=list(
    'rootDir' = "/home/kimlab/pancreas.plasma.ev.long.RNA",
    'findFilesFunc' = findPancreasPlasmaEV_longRNAFiles, 
    'getColDFFunc' = getPancreasPlasmaEV_longRNAColDataDF,        
    'tx2MappingFile'='gencode.v35.ucsc.rmsk.tx.to.gene.csv',
    'tx2MappingDir' = file.path('/home/kimlab/genomes.annotations/gencode.35'),
    'outputRoot' = "/home/aedavids/extraCellularRNA/data/R/output",
    #'outputRoot' = "/home/kimlab/aedavids/extraCellularRNA/data/R/output",
    'mapToIdName'='biotype',
    'normalizedFileName'='pancreas.plasma.ev.long.RNA.normalized.deseq',
    'transcriptMapperDFFunc'=getTETBioTypeMapper,
    'ignoreAfterBar'=FALSE,
    'txOut'=FALSE,
    'designStr' = "~ diseaseState"
  )
)

missingParmeters <- function( warningLabel, parameterList) {
  #cat( "AEDWIP ", names(parameterList))
  error = FALSE
  keyList = list('rootDir', 'findFilesFunc', 'getColDFFunc', 'tx2MappingFile', 
                 'tx2MappingDir', 'mapToIdName', 'normalizedFileName', 
                 'outputRoot', 'transcriptMapperDFFunc',
                 'ignoreAfterBar', 'txOut', 'designStr')
  
  for (key in keyList) {
    if (! exists(key, parameterList) ) {
      error = TRUE
      # warning( paste(warningLabel, key, " is missing") )
      warning( paste("\n", warningLabel, key," is missing")  )
    }
  }
  
  # find unexpected parameters. probably spelling error
  badKeys = setdiff( names(parameterList), keyList )
  if ( length(badKeys) != 0 ) {
    warning(paste("\n", warningLabel, "unknown keys", badKeys))
  }
  
  return( error )
}

# test driver
testSet = 'debug'
isMissing = missingParmeters(testSet, debugParametersList[[testSet]])
if (isMissing) {
  warning('\nERROR')
}
