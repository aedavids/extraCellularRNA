
#
# TE.R
# functions for parsing and working with gencode.v35.ucsc.rmsk.tx.to.gene.csv
#
# aedavids@ucsc.edu
#

library('dplyr')
library('readr')

########################################################################
getTETxMapper <- function(tx2MappingDir, tx2MappingFile) {
  # function passes in parameter list
  # enable lazy evaluation
  tx2MappingFilePath <- file.path(tx2MappingDir, tx2MappingFile)
  ret <- TE_TranscriptMapper(tx2MappingFilePath)
  
  return (ret)
}

########################################################################
getTETBioTypeMapper <- function(tx2MappingDir, tx2MappingFile) {
  # function passes in parameter list
  # enable lazy evaluation
  tx2MappingFilePath <- file.path(tx2MappingDir, tx2MappingFile)
  retDF <- TE_bioTypeTranscriptMapper(tx2MappingFilePath) 
  
  return (retDF)
}

########################################################################
getTEGeneMapper <- function(tx2MappingDir, tx2MappingFile) {
  # function passes in parameter list
  # enable lazy evaluation
  
  # warning(sprintf("getTEGeneMapper tx2MappingDir:%s tx2MappingFile:%s \n", 
  #                 tx2MappingDir, tx2MappingFile), immediate.=TRUE)
  
  tx2MappingFilePath <- file.path(tx2MappingDir, tx2MappingFile)
  
  # warning(sprintf("getTEGeneMapper() tx2MappingFilePath:%s \n", 
  #                 tx2MappingFilePath), immediate.=TRUE)
  # 
  # warning(sprintf( "getTEGeneMapper() file.exists() %s \n",
  #                  as.character( file.exists(tx2MappingFilePath))), 
  #         immediate.=TRUE)
  
  retDF <- TE_geneTranscriptMapper(tx2MappingFilePath)
  
  return (retDF)
}

########################################################################
TE_TranscriptMapper <- function(tx2MappingFilePath) {
  # returns a 2 column data frame that can be used with tximport().
  # purpose is to enable DSEQ count normalization at the transcript level.
  # both columns will be the same.
  #
  # A tibble: 6 x 2
  # tx                                                                                        biotype 
  # <chr>       
  #   5 ENST00000473358.1|ENSG00000243485.5|OTTHUMG00000000959.2|OTTHUMT00000002840.1|MIR1302-… MIR1302-2…
  #   1 hg38_rmsk_L1PA2_range=chr22_KI270739v1_random:41877-42470_5'pad=0_3'pad=0_strand=+_repea… L1PA2   
  #   2 hg38_rmsk_Tigger5b_range=chr22_KI270739v1_random:42472-42520_5'pad=0_3'pad=0_strand=-_re… Tigger5b
  #
  # motivation:
  # we can not use generic mapping files like
  # kimlab/genomes.annotations/gencode.35/gencode.v35.tx.to.gene.csv
  # each row is of the following format
  # ENST00000488147.1|ENSG00000227232.5|OTTHUMG00000000958.1|OTTHUMT00000002839.1|WASH7P-201|WASH7P|1351|unprocessed_pseudogene|,WASH7P
  # a row is a single string of ids seperated by "|".
  # ref "
  # https://github.com/rreggiar/welcome-to-the-kim-lab/wiki/Salmon-output-data-dictionary#understanding-the-name-column
  #
  # 'gencode.v35.ucsc.rmsk.tx.to.gene.csv' does not follow this format
  # there are two columns
  # the second column is the gene name
  # the first is to be used as the transcript id
  # The transcript ids at the top of the file follow the above format
  #
  # the transcripts at the bottom represent TE's and are of the form
  # hg38_rmsk_ALR/Alpha_range=chr22_KI270739v1_random:777-13807_5'pad=0_3'pad=0_strand=+_repeatMasking=none,ALR/Alpha
  # hg38_rmsk_ALR/Alpha_range=chr22_KI270739v1_random:14128-35940_5'pad=0_3'pad=0_strand=+_repeatMasking=none,ALR/Alpha
  # hg38_rmsk_L1PA2_range=chr22_KI270739v1_random:35954-37185_5'pad=0_3'pad=0_strand=+_repeatMasking=none,L1PA2
  #
  # they do not encode multiple ids. 
  #
  
  DF <- read_csv(tx2MappingFilePath, col_names=c('tx', 'gene'))
  c1 <- DF[,1]
  retDF <- cbind(c1, c1)
  
  return( retDF )
}

########################################################################
TE_geneTranscriptMapper <- function(tx2MappingFilePath) {
  # returns a 2 column data frame that can be used with tximport().
  
  DF <- read_csv(tx2MappingFilePath, col_names=c('tx', 'gene'))
  return( DF )
}

########################################################################
TE_bioTypeTranscriptMapper <- function(tx2MappingFilePath) {
  # returns a 2 column data frame that can be used with tximport().
 
  # AEDWIP TODO root only works in Andy's docker image
  root <- '/home/kimlab'
  te.family.clade.csvPath <- file.path( root,
        'exoRNA-biomarkers-panc/output.data/reference.data/te.family.clade.csv')
  
  # AEDWIP TODO create a second version of this function that uses 'cladeBiotype'
  retDF <- getRMASKMapDF(tx2MappingFilePath, te.family.clade.csvPath, 
                         teBiotype='familyBiotype')
  
  # The trascript id will be something like
  # DDX11L1-202 for standard format salmon quant.sf transcripts
  # ALR/Alpha_chr22_5'pad=0 for  TE formated transcripts 
  # DESeq2 can not process this format because the tx id does not string
  # match the quant.sf entries
  # we need the entire first column of the tx2GeneDF
  retDF <- retDF[,c('txLong', 'biotype')]
  
  return( retDF )
}

############################################################
getRMASKMapDF <- function(tx2MappingFilePath, te.family.clade.csvPath, 
                          teBiotype='familyBiotype') {
  # AEDWIP TODO
  # this is really slow save/load from cache see testTEParsing.Rmd
  #
  # see testTEParsing.Rmd
  #
  # arguments
  #   tx2MappingFilePath
  #     kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv
  #
  #   te.family.clade.csvPath
  #     kimlab/exoRNA-biomarkers-panc/output.data/reference.data/te.family.clade.csv
  #
  #   teBiotype
  #     either 'familyBiotype' or 'cladeBiotype'
  #     default is 'familyBiotype'
  #     TE's do not conform to conventional biotypes. clade or family is use
  #     familyBiotype, more general (4 total types)
  #     cladeBiotype clade (dozens of types)
  
  
  if (!(teBiotype == "cladeBiotype" | teBiotype == "familyBiotype")) {
    fmt <- "unknown teBiotype:%s must be 'cladeBiotype' or 'familyBiotype'"
    errMsg <- sprintf(fmt, as.character(te))
    stop(errMsg)
  }
  
  # readr read_table if faster than built in
  # mapDF <- read_table(tx2MappingFilePath)
  mapDF <- read_csv( tx2MappingFilePath, col_names=c( 'id', 'gene') )
  #colnames(mapDF) <- c('gene')
  
  # readr read_csv if faster than built in
  te.family.cladeDF <- read_csv(te.family.clade.csvPath, 
                                col_names=c('gene', 'family', 'clade'),
                                col_type=c(col_character(), col_character(), col_character())
  ) %>%
    distinct()
  
  stdMapDF <- parseStandardFormat_(mapDF)
  
  retColNames <- c('txLong', 'gene', 'biotype')
  stdMapDF <- stdMapDF[,retColNames]
  
  if (teBiotype == 'familyBiotype') {
    selectCols <- c('txLong', 'gene', 'familyBiotype')
  } else {
    selectCols <- c('txLong', 'gene', 'cladeBiotype')
  }
  teMapDF <- parseTEFormat_(mapDF, te.family.cladeDF)
  teMapDF <- teMapDF[, selectCols]
  names(teMapDF) <- retColNames
  
  retDF <- rbind(stdMapDF, teMapDF )
  
  # # add the first column as the 'long tx name'
  # # DESeq2 can not handle our TE transcript id, it does not follow
  # # expected format
  # txLong <- mapDF[, c('id')]
  # retDF <- cbind(retDF, txLong)

  return( retDF )  
}

########################################################################
parseStandardFormat_ <- function(mapDF, debug=FALSE) {
  # returns a data frame with columns tx, gene, len, biotype
  
  # ref:
  # mappingFile <- 'gencode.v35.ucsc.rmsk.tx.to.gene.csv'
  # indexFilePath <- file.path( '/home/kimlab/genomes.annotations/gencode.35', mappingFile)
  
  # this is slow, try using read_csv, it returns tibble instead of a dataframe
  # mapDF <- read.table(indexFilePath, header=FALSE, col.names=c('gene'))
  # AEDWIP TODO: 2 col data frame , id, gene
  
  # parse the top of the index file
  # kimlab/genomes.annotations/gencode.35/gencode.v35.tx.to.gene.csv
  # each row is of the following format
  # ENST00000488147.1|ENSG00000227232.5|OTTHUMG00000000958.1|OTTHUMT00000002839.1|WASH7P-201|WASH7P|1351|unprocessed_pseudogene|,WASH7P
  # a row is a single string of ids seperated by "|".
  # ref "
  # https://github.com/rreggiar/welcome-to-the-kim-lab/wiki/Salmon-output-data-dictionary#understanding-the-name-column
  #  
  # Roman wrote this code
  # gen <- read_csv(paste0(output.data.dir, 
  #                      path, '/',                                
  #                     paste0(formula, '/'), 'tx.',ref,'.',path,'.', 
  #                           first,'.v.',second,'.de-seq.counts.csv')
  #                     ) %>% 
  # rename(X1 = 'gene') %>% 
  # filter_at(vars(contains('.')), any_vars(. >= 10)) %>% 
  # filter(grepl('^ENST', gene)) %>% 
  # separate(gene, sep = '\\|', into = c(NA,NA,NA,NA,'tx','gene','len','biotype'))
  
  
  #warning("AEDWIP begin parseTop_()\n")
  
  #warning("AEDWIP DEBUG add filter()\n")
  #headMapDF <- mapDF %>% 
    #rename('V1'='gene') %>% 
    
    # TODO what does filter_at do? what is example in and out?
    # https://dplyr.tidyverse.org/reference/filter_all.html
    # https://dplyr.tidyverse.org/reference/vars.html
    #filter_at(vars(contains('.')), any_vars(. >= 10)) %>% 
    
    # filter() function is used to subset a data frame, retaining all rows that satisfy your conditions
    # https://dplyr.tidyverse.org/reference/filter.html    
    # grepl returns a logical vector (match or not for each element of x).
    #filter(grepl('^ENST', gene)) %>% 
    ff <- filter(mapDF, grepl('^ENST', id))
    if (debug) {
      print("after filter")
      print(ff)
    }
    
    # keep a copy of the original id
    # DESeq2 can not handle our TE short tx ids
    mm <- mutate(ff, 'txLong' = id) 
    if (debug) {
      print("after mutate")
      print(mm)
    }
    
    # ignore warning
    # Expected 8 pieces. Additional pieces discarded in 6 rows [1, 2, 3, 4, 5, 6].
    ss <- tidyr::separate(mm, id, sep="\\|", into=c(NA, NA, NA, NA, 'tx', 'gene', 'len', 'biotype'))
    if (debug) {
      print("after separate")
      print(ss)
    }
    
    # select( tx, gene, biotype, txLong)
    
    retDF <- ss
  
  #warning("\nAEDWIP END parseTop_()\n")
  
  return( retDF)
}

########################################################################
parseTEFormat_ <- function(tx2GeneDF, te.family.cladesDF, debug=FALSE) {
  #
  # parse Roman's custom TE format, see examples in comment bellow
  #
  # arguments
  #   tx2GeneDF,
  #     2 columns, 'id', and 'gene'
  #         gene family clade
  #       1   L1M2     L1  LINE
  #       2 AluSx3    Alu  SINE
  #       3  AluSx    Alu  SINE
  #       4  AluJr    Alu  SINE
  #
  # debug
  #   boolean if true prints out results after each pipeline stp
    
  # filter_at() filters rows in select columns. The syntax is pretty garbage, 
  # but the first argument is vars() that you provide your filtering criteria 
  # too (contains(‘string’), is.numeric(), etc) the second is the criteria you 
  # want to filter with — any_vars(.  > 10). 
  # Error: `.predicate` has no matching columns.
  # fa <- filter_at(tx2GeneDF, vars(contains('.')), any_vars(. >= 10))
  # if (debug) {
  #   print("after filter_at()")
  #   print(fa)
  # }
    
  # remove rows that start with ENST
  #    filter(!grepl('^ENST', gene)) %>% 
  f <- filter(tx2GeneDF, !grepl('^ENST', id)) 
  if (debug) {
    print("after filter")
    print(f)
    # [1] "hg38_rmsk_L1PA2_range=chr22_KI270739v1_random:41877-42470_5'pad=0_3'pad=0_strand=+_repeatMasking=none"    
    # [2] "hg38_rmsk_Tigger5b_range=chr22_KI270739v1_random:42472-42520_5'pad=0_3'pad=0_strand=-_repeatMasking=none" 
    # [3] "hg38_rmsk_ALR/Alpha_range=chr22_KI270739v1_random:42531-48656_5'pad=0_3'pad=0_strand=-_repeatMasking=none"
    # [4] "hg38_rmsk_L1PA3_range=chr22_KI270739v1_random:48657-49096_5'pad=0_3'pad=0_strand=-_repeatMasking=none"    
    # [5] "hg38_rmsk_L1PA3_range=chr22_KI270739v1_random:49092-50070_5'pad=0_3'pad=0_strand=+_repeatMasking=none"    
    # [6] "hg38_rmsk_ALR/Alpha_range=chr22_KI270739v1_random:50071-73985_5'pad=0_3'pad=0_strand=-_repeatMasking=none"
  }
  
  # keep a copy of the original id
  # DESeq2 can not handle our TE short tx ids
  mi <- mutate(f, 'txLong' = id)
    
  # ignore warning
  # Expected 8 pieces. Additional pieces discarded in 6 rows [1, 2, 3, 4, 5, 6].
  # https://www.rdocumentation.org/packages/tidyr/versions/0.8.3/topics/separate
  #   separate(gene, sep = '_', into = c(NA, NA, 'gene', 'position', NA, NA, 
  #                                     'strand', NA)) %>% 
  s <- tidyr::separate(mi, id, sep = '_', into=c(NA, NA, 'gene', 'position', NA, NA,
                                           'strand', NA))
  if (debug) {
    print("after separate")
    print(s)
    #     gene      position    strand 
    #     <chr>     <chr>       <chr>  
    #   1 L1PA2     range=chr22 5'pad=0
    #   2 Tigger5b  range=chr22 5'pad=0  
  }

  mu <- mutate(s, position = sub('range=', '', position),
           strand = sub('strand=', '', strand))
  if (debug) {
    print("after mutate")
    print(mu)
    # A tibble: 6 x 3
    # gene      position strand 
    # <chr>     <chr>    <chr>  
    #   1 L1PA2     chr22    5'pad=0
    # 2 Tigger5b  chr22    5'pad=0  
  }

  # AEDWIP merge does join. are we droping rows?
  me <- merge(mu, te.family.cladesDF %>% distinct(), by = 'gene') 
  if (debug) {
    print("after merge")
    print(me)
    #   gene         position  strand     family clade
    # 1 ALR/Alpha    chr22 5'pad=0        centr Satellite
    # 2 ALR/Alpha    chr22 5'pad=0        centr Satellite 
  }
    
  mu2 <- mutate(me, tmp.gene = gene) 
  if (debug) {
    print("after second mutate")
    print(mu2)
    #    gene        position  strand       family     clade  tmp.gene
    # 1 ALR/Alpha    chr22 5'pad=0        centr Satellite ALR/Alpha
    # 2 ALR/Alpha    chr22 5'pad=0        centr Satellite ALR/Alpha 
    }

  # unite() Convenience function to paste together multiple columns into one.
  # https://tidyr.tidyverse.org/reference/unite.html
  u <- tidyr::unite(mu2, 'tx', c(gene, position, strand), sep = '_') 
  if (debug) {
    print("after unite")
    print(u)
    #     tx                             family     clade  tmp.gene
    #   1 ALR/Alpha_chr22_5'pad=0        centr Satellite ALR/Alpha
    #   2 ALR/Alpha_chr22_5'pad=0        centr Satellite ALR/Alpha
    #   3     L1PA2_chr22_5'pad=0           L1      LINE     L1PA2
  }
  
  # Roman " TE’s have a dual-encoded biotype —> Family (len) and Clade 
  # (biotype). You can drop len and just use clade as biotype. Family is more 
  # general (4 total) while clade has a bit more resolution (dozens)"
  #
  ret <- select(u, tx, 'gene' = tmp.gene, 'cladeBiotype' = clade, 
                'familyBiotype' = family, txLong, everything())
  
  if (debug) {
    print("after select")
    print(ret)
    #   tx                      gene       cladeBiotype      familyBiotype
    # 1 ALR/Alpha_chr22_5'pad=0 ALR/Alpha Satellite        centr
    # 2 ALR/Alpha_chr22_5'pad=0 ALR/Alpha Satellite        centr
    # 3     L1PA2_chr22_5'pad=0     L1PA2      LINE           L1  
  }
  return(ret)
}
