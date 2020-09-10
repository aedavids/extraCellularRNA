

#
# function to find salmon quant.sf files
#
########################################################################
findBulkDay5CntrlFiles <- function (bulkDirPath) {
  #bulkDirPath
  
  f <- list.files( bulkDirPath )
  cntrlIdxList = grep('ctrl', f)
  ctrlDirs = f[ cntrlIdxList ]
  
  bulkDay5CntrlSalmonFiles <- paste( bulkDirPath, ctrlDirs, 
                                     'gencode.salmon.out/quant.sf', 
                                     sep="/" )
  return( bulkDay5CntrlSalmonFiles )
}

########################################################################
findBulkDay5KrasFiles <- function (bulkDirPath) {
  f <- list.files( bulkDirPath )
  kraslIdxList = grep('kras', f)
  krasDirs = f[ kraslIdxList ]
  #krasDirs
  
  bulkDay5KrasSalmonFiles <- paste( bulkDirPath, krasDirs, 
                                    'gencode.salmon.out/quant.sf',
                                    sep="/" )
  
  # [4] "/home/kimlab/kras.ipsc/bulk.data/day.5/kras.alu.edit.out/gencode.salmon.out/quant.sf"
  bulkDay5KrasSalmonFiles <- bulkDay5KrasSalmonFiles[1:3]
  return( bulkDay5KrasSalmonFiles )
}

########################################################################
findExoDay5CtrlFiles <- function( exoDirPath ) {
  #list.dirs( exoDirPath )
  f = list.files( exoDirPath )
  
  ctrlIdxList = grep('ctrl', f)
  cntrlDirs = f[ ctrlIdxList ]
  cntrlDirs
  
  exoDay5CntrlSalmonFiles <- paste( exoDirPath, cntrlDirs, 
                                    'gencode.salmon.out/quant.sf', 
                                    sep="/" )
  exoDay5CntrlSalmonFiles
}

########################################################################
findExoDay5KrasFiles <- function( exoDirPath ) {
  f = list.files( exoDirPath )
  kraslIdxList = grep('kras', f)
  krasDirs = f[ kraslIdxList ]
  krasDirs
  
  exoDay5KrasSalmonFiles <- paste( bulkDirPath, krasDirs, 
                                   'gencode.salmon.out/quant.sf', 
                                   sep="/" )
  return( exoDay5KrasSalmonFiles )
}

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
checkIndexIsSame <- function( cmdJSONList, refIdxFilePath) {
  # returns true if all the index files are the same
  
  idxIsSame = lapply(cmdJSONList, function(filePath){ filePath == refIdxFilePath })
  idxIsSame <- unlist( idxIsSame )
  return( all(idxIsSame) )
}
