## Part 2. Search for G4 loci that correspond to each TF. 


library(tidyr)
library(dplyr)
library(readr)
library(parallel)
library(openxlsx)
library(GenomicRanges)
library(ArchR)
library(EpiTrace)
library(BSgenome.Mmusculus.UCSC.mm10)

dir.create('/gpfs/output/ECS_Research/data/G4_Aging/DNAm/TFBS')
setwd('/gpfs/output/ECS_Research/data/G4_Aging/DNAm/TFBS')

mm10_g4_positive_all_sites <- read_tsv('/gpfs/output/ECS_Research/data/G4_Aging/bed/ref_all/mm10_g4_positive_all_sites.bed',col_names=c('chr','start','end')) %>% makeGRangesFromDataFrame()
mm10_g4_negative_all_sites <- read_tsv('/gpfs/output/ECS_Research/data/G4_Aging/bed/ref_all/mm10_g4_negative_all_sites.bed',col_names=c('chr','start','end')) %>% makeGRangesFromDataFrame()
mm10_g4_positive_aging_sites <- read_tsv('/gpfs/output/ECS_Research/data/G4_Aging/bed/ref_all/mm10_g4_positive_aging_open_sites.bed',col_names=c('chr','start','end')) %>% makeGRangesFromDataFrame()
mm10_g4_negative_aging_sites <- read_tsv('/gpfs/output/ECS_Research/data/G4_Aging/bed/ref_all/mm10_g4_negative_aging_open_sites.bed',col_names=c('chr','start','end')) %>% makeGRangesFromDataFrame()

chrom_lengths <- seqlengths(BSgenome.Mmusculus.UCSC.mm10) - 8

list_of_interested_regions <- list(
  'pGQSpos' = mm10_g4_positive_all_sites,
  'pGQSneg' = mm10_g4_negative_all_sites,
  'pGQSpos_aging' = mm10_g4_positive_aging_sites,
  'pGQSneg_aging' = mm10_g4_negative_aging_sites
)

lapply(list_of_interested_regions,function(gr){
  gr <- gr[gr@seqnames %in% names(chrom_lengths)]
  subsetByOverlaps(gr, GRanges(seqnames = names(chrom_lengths), 
                               ranges = IRanges(start = 1, end = chrom_lengths))) 
}) -> clean_list_of_interested_regions
names(clean_list_of_interested_regions) <- names(list_of_interested_regions)

motifs <- readRDS('/gpfs/output/ECS_Research/data/G4_Aging/DNAm/TFBS/JASPAR2022.motifs.mm10.rds') 
motifSummary <- readRDS('//gpfs/output/ECS_Research/data/G4_Aging/DNAm/TFBS/JASPAR2022.motifSummary.mm10.rds') 

mclapply(clean_list_of_interested_regions,mc.cores = 30,function(peakSet){
  cutOff=5e-05
  width=7
  BSgenome = 'BSgenome.Mmusculus.UCSC.mm10'
  motifPositions <- motifmatchr::matchMotifs(pwms = motifs, 
                                             subject = peakSet, genome = BSgenome, out = "positions", 
                                             p.cutoff = cutOff, w = width)
  motifPositions
}) -> calculated_motif_positions



lapply(names(calculated_motif_positions),function(x){
  clean_list_of_interested_regions[[x]] -> peak_regions
  calculated_motif_positions[[x]] -> motif_regions_list 
  mclapply(names(motif_regions_list),mc.cores = 40,function(y){
    peak_regions[findOverlaps(motif_regions_list[[y]],peak_regions)@to%>%unique(),]
  }) -> return_peak_regions_for_each_motif
  names(return_peak_regions_for_each_motif) <- names(motif_regions_list)
  return_peak_regions_for_each_motif 
}) -> calculated_motif_overlapping_peak_regions
names(calculated_motif_overlapping_peak_regions) <- names(calculated_motif_positions)


saveRDS(calculated_motif_positions,file='TFBS.positions.for.each.peakSet.mm10.rds')
saveRDS(clean_list_of_interested_regions,file='peakSet.list.for.TFBS.DNAm.analysis.mm10.rds')
saveRDS(calculated_motif_overlapping_peak_regions,file='TFBS.overlapping.peak.positions.for.each.peakSet.mm10.rds')
