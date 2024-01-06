## Note: the Part 1-4 of human G4AMP analysis is similar to mouse, except that:
### a) exchanging the mouse genome (mm10) to human (hg19) genome in "Part 1", and 
### b) using human G4 peak set instead of mouse G4 peak set in "Part 2". 

# Detail analysis of TF-G4 DNAm 

library(dplyr)
library(tidyr)
library(reshape2)
library(openxlsx)
library(ArchR)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(Seurat)
library(readr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# saveRDS(G4_region_TFBS_beta_df_summarised_for_t3_new %>% dplyr::select(library_name,real_age,Group,Mutated_gene) %>% unique(),file='/Users/wing/Desktop/Eulerian/文件/内部研发项目/G4_aging/data/G4AMP_progeroid/meta_human.rds')

meta_human <- readRDS('/Users/wing/Desktop/Eulerian/文件/内部研发项目/G4_aging/data/G4AMP_progeroid/meta_human.rds')
human_rdslist <- list.files(path='/Users/wing/Desktop/Eulerian/文件/内部研发项目/G4_aging/data/G4AMP_progeroid/rds',full.names = T,pattern='*rds')

library_id_ctrl <- meta_human$library_name[meta_human$Group %in% 'CTRL']
library_id_progeroid <- meta_human$library_name[meta_human$Group %in% 'Progeroid']

ctrl_rdslist <- human_rdslist[grepl(paste0(library_id_ctrl,collapse = '|'),human_rdslist)]
progeroid_rdslist <- human_rdslist[grepl(paste0(library_id_progeroid,collapse = '|'),human_rdslist)]

lapply(ctrl_rdslist,function(x){
  readRDS(x) -> y 
  y$tfbs_result -> z
  z$library_name <- gsub('\\..+','',basename(x))
  z
}) %>% bind_rows() -> G4_region_TFBS_beta_df_for_ctrl

G4_region_TFBS_beta_df_for_ctrl <- left_join(G4_region_TFBS_beta_df_for_ctrl,meta_human)

lapply(progeroid_rdslist,function(x){
  readRDS(x) -> y 
  y$tfbs_result -> z
  z$library_name <- gsub('\\..+','',basename(x))
  z
}) %>% bind_rows() -> G4_region_TFBS_beta_df_for_progeroid

G4_region_TFBS_beta_df_for_progeroid <- left_join(G4_region_TFBS_beta_df_for_progeroid,meta_human)

G4_region_TFBS_beta_df_for_ctrl %>% dplyr::group_by(TF) %>% rstatix::cor_test(vars = c('beta','real_age')) -> human_G4_vs_age_cor_res
human_G4_vs_age_cor_res$padj <- p.adjust(human_G4_vs_age_cor_res$p)

G4_region_TFBS_beta_df_for_ctrl %>% dplyr::filter(TF %in% (human_G4_vs_age_cor_res %>% dplyr::filter(padj<0.001 & cor < -0.5) %>% dplyr::select(TF) %>% unlist())) %>% 
  dplyr::group_by(library_name,real_age,Group,Mutated_gene) %>% 
  dplyr::summarise(G4_TF_beta=sum(coverage*beta)/sum(coverage)) -> G4_region_TFBS_beta_df_summarised_for_CTRL

G4_region_TFBS_beta_df_for_progeroid %>% dplyr::filter(TF %in% (human_G4_vs_age_cor_res %>% dplyr::filter(padj<0.001 & cor < -0.5) %>% dplyr::select(TF) %>% unlist())) %>% 
  dplyr::group_by(library_name,real_age,Group,Mutated_gene) %>% 
  dplyr::summarise(G4_TF_beta=sum(coverage*beta)/sum(coverage)) -> G4_region_TFBS_beta_df_summarised_for_progeroid

G4_region_TFBS_beta_df_summarised_combined <- rbind(G4_region_TFBS_beta_df_summarised_for_CTRL,G4_region_TFBS_beta_df_summarised_for_progeroid)


ggplot(G4_region_TFBS_beta_df_summarised_combined,aes(x=real_age,y=G4_TF_beta)) + scale_color_manual(values=c('gray','orange')) + geom_point(pch=21,aes(fill=Mutated_gene),size=6,color='black')  + geom_smooth(method='lm',se=F,aes(group=Group,linetype=Group,color=Group)) + scale_fill_manual(values=c('pink','beige','red')) + theme_classic() + theme(text=element_blank(),legend.position = 'none') + scale_alpha_manual(values=c(0.01,0.5))  + ylab('\nG4 region DNAm\n') + scale_linetype_manual(values=c('solid','dashed')) + ggpubr::stat_cor(data=G4_region_TFBS_beta_df_summarised_combined %>% dplyr::filter(library_name %ni% '1902756S1G4'), mapping=aes(group=Group),label.y.npc = 1,label.x.npc = 0.5)
