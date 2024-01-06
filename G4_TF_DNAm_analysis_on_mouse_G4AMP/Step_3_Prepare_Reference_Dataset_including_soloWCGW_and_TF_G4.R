## Part 3. Prepare reference data list (including solo_WCGW CpG and TF-bound G4)

library(tidyr)
library(dplyr)
library(readr)
library(parallel)
library(openxlsx)
library(GenomicRanges)
library(ArchR)
library(patchwork)


# Clock loci
read_tsv('/gpfs/output/ECS_Research/data/G4_Aging/DNAm/non_TFBS/solo_WCGW_inCommonPMDs_mm10.bed.gz',col_names=c('chr','start','end')) %>% makeGRangesFromDataFrame() -> ref_solo_wcgw_gr

## TFBS DNAm

calculated_motif_overlapping_peak_regions <- readRDS('/gpfs/output/ECS_Research/data/G4_Aging/DNAm/TFBS/TFBS.overlapping.peak.positions.for.each.peakSet.mm10.rds')
TF_negative_union_list <- readRDS('/gpfs/output/G4AMP/bin/TF_negative_union_selected.20230522.rds')


DNAm_analysis_reference_data = list(
  clock_reference = list(
    'Clock_Solo_WCGW_gr' = ref_solo_wcgw_gr
  ),
  TFBS_reference = list(
    'pGQSpos' = calculated_motif_overlapping_peak_regions[['pGQSpos']],
    'pGQSneg' = calculated_motif_overlapping_peak_regions[['pGQSneg']],
    'pGQSpos_aging' = calculated_motif_overlapping_peak_regions[['pGQSpos_aging']],
    'pGQSneg_aging' = calculated_motif_overlapping_peak_regions[['pGQSneg_aging']],
    'TF_tumor_specific_negative_selected' = TF_negative_union_list
  )
)

saveRDS(DNAm_analysis_reference_data,file='/gpfs/output/ECS_Research/data/G4_Aging/DNAm/DNAm_analysis_reference_data.for.mm10.rds')