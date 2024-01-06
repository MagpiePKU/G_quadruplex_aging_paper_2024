## Part 1. Extract TFBS Motifs


library(tidyr)
library(dplyr)
library(readr)
library(parallel)
library(openxlsx)
library(GenomicRanges)
library(ArchR)
library(patchwork)
library(EpiTrace)
library(BSgenome.Mmusculus.UCSC.mm10)

collection = 'CORE'
genome_ref = 'mm10'
species = 'Mus musculus'


summarizeJASPARMotifs <- function(motifs = NULL){ # copied from ArchR. 
  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)
  motifDF <- lapply(seq_along(motifs), function(x){
    data.frame(
      row.names = motifNames[x],
      name = motifs[[x]]@name[[1]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      symbol = ifelse(!is.null(motifs[[x]]@tags$symbol[1]), motifs[[x]]@tags$symbol[1], NA) ,
      family = ifelse(!is.null(motifs[[x]]@tags$family[1]), motifs[[x]]@tags$family[1], NA),
      alias = ifelse(!is.null(motifs[[x]]@tags$alias[1]), motifs[[x]]@tags$alias[1], NA),
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame
  names(motifs) <- motifNames
  out <- list(motifs = motifs, motifSummary = motifDF)
  return(out)
}

args <- list(species = species, collection = collection)
motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, args) # worked. now get 633 motifs MA0030.1 etc. 
motifs <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, args)
obj <- summarizeJASPARMotifs(motifs)
motifs <- obj$motifs
motifSummary <- obj$motifSummary

dir.create('/Users/wing/Desktop/Eulerian/文件/内部研发项目/G4_aging/data/G4_TFBS_analysis')
saveRDS(motifSummary,file='/Users/wing/Desktop/Eulerian/文件/内部研发项目/G4_aging/data/G4_TFBS_analysis/JASPAR2022.motifSummary.mm10.rds')
saveRDS(motifs,file='/Users/wing/Desktop/Eulerian/文件/内部研发项目/G4_aging/data/G4_TFBS_analysis/JASPAR2022.motifs.mm10.rds')
