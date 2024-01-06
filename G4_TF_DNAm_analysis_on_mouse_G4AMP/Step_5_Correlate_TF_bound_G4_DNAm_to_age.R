## Part 5. Correlate TF-G4 DNAm to age

library(tidyr)
library(dplyr)
library(readr)
library(parallel)
library(openxlsx)
library(GenomicRanges)
library(ArchR)
library(patchwork)


rdslist <- list.files(path='/Users/wing/Desktop/Eulerian/文件/内部研发项目/G4_aging/data/G4AMP_mouse/20240101_mouse_G4AMP_aging/rds',pattern='*rds$',full.names = T)

lapply(rdslist,function(x){
  readRDS(x) 
}) -> rds_res_list

names(rds_res_list) <- gsub('G4.DNAm.simple.output.rds','',basename(rdslist)) 

lapply(names(rds_res_list),function(x){
  rds_res_list[[x]]$tfbs_result -> y
  y$source <- x
  y
}) %>% bind_rows() -> result_DNAm

read.xlsx('/Users/wing/Desktop/Eulerian/文件/内部研发项目/G4_aging/data/G4AMP_mouse/Mouse_meta/Mouse Sample List.xlsx',sheet=1) -> meta_mouse
result_DNAm$Family.ID <- gsub('-.+','',result_DNAm$source)
result_DNAm <- left_join(result_DNAm,meta_mouse)
result_DNAm$age <- result_DNAm$`Age.(day)`

result_DNAm %>% dplyr::filter(regionSet %in% 'pGQSpos') %>% dplyr::group_by(TF) %>% rstatix::cor_test(vars=c('beta','age')) -> age_x_beta_pGQSpos

result_DNAm %>% dplyr::filter(TF %in% (age_x_beta_pGQSpos %>% dplyr::filter(padj<0.01 & cor < 0.5) %>% dplyr::select(TF) %>% unlist()) & regionSet %in% 'pGQSpos') %>% 
  dplyr::group_by(source,age) %>% 
  dplyr::summarise(beta=sum(coverage*beta)/sum(coverage)) -> result_DNAm_aging_TF

age_x_beta_pGQSpos$padj <- p.adjust(age_x_beta_pGQSpos$p)

ggplot(age_x_beta_pGQSpos,aes(x=cor)) + geom_vline(xintercept = 0,linetype='dashed')  + geom_histogram() + theme_classic() + theme(text=element_text(size=20)) + xlab('\nCorrelation\nage x TF-positive peak DNAm') + ylab('\n\n\nTF Count\n') -> p1
ggplot(result_DNAm %>% dplyr::filter(TF %in% 'Stat2_80' & regionSet %in% 'pGQSpos'),aes(x=age,y=beta)) + geom_point() + scale_x_log10() + ggpubr::stat_cor(size=6,label.y.npc = 0.1) + geom_smooth(method='lm') + theme_classic() + theme(text=element_text(size=20)) + ylab('\n\n\nStat2-positive G4 DNAm\n') + xlab('\nlog Sample age (days)\n') -> p2
ggplot(result_DNAm %>% dplyr::filter(TF %in% 'Irf1_104' & regionSet %in% 'pGQSpos'),aes(x=age,y=beta)) + geom_point() + scale_x_log10() + ggpubr::stat_cor(size=6,label.y.npc = 0.1) + geom_smooth(method='lm') + theme_classic() + theme(text=element_text(size=20)) + ylab('\n\n\nIrf1-positive G4 DNAm\n') + xlab('\nlog Sample age (days)\n') -> p3
ggplot(result_DNAm_aging_TF,aes(x=age,y=beta)) + geom_point() + scale_x_log10() + ggpubr::stat_cor(size=6,label.y.npc = 0.1) + geom_smooth(method='lm') + theme_classic() + theme(text=element_text(size=20)) + ylab('\nAge-associated TF G4 DNAm\n') + xlab('\nlog Sample age (days)\n') -> p4

wrap_plots(list(p1,p2,p3,p4),nrow=1)