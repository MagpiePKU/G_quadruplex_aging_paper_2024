## Part 4. Calculate solo-WCGW DNAm and TF-bound G4 DNAm for each sample 

Calculate_Clock_DNAm_for_G4AMP <- function(input_methylKit_df,input_methylKit_gr,DNAm_analysis_reference_data_list=DNAm_analysis_reference_data){
  
  y = input_methylKit_df
  y_gr = input_methylKit_gr
  
  # solo
  y_overlap <- findOverlaps(y_gr,DNAm_analysis_reference_data_list$clock_reference$Clock_Solo_WCGW_gr)
  y_comb <- y[unique(y_overlap@from),]
  solo_out=sum(y_comb$freqC*y_comb$coverage)/sum(y_comb$coverage)
  solo_cov = sum(y_comb$coverage)
  solo_n = nrow(y_comb)
  
  # total cpg coverage QC
  y$coverage %>% sum -> total_cpg_cov
  y$coverage %>% length -> total_cpg_n
  sum(y$freqC*y$coverage)/sum(y$coverage) -> total_cpg_beta
  
  return_clock_dnam_result_list = list(
    solo_out=solo_out,
    solo_cov=solo_cov,
    solo_n = solo_n,
    total_cpg_beta=total_cpg_beta,
    total_cpg_cov=total_cpg_cov,
    total_cpg_n=total_cpg_n
  )
  
  return(return_clock_dnam_result_list)
  
}

Calculate_TFBS_DNAm_for_G4AMP <- function(input_methylKit_df,input_methylKit_gr,DNAm_analysis_reference_data_list=DNAm_analysis_reference_data){
  
  y = input_methylKit_df
  y_gr = input_methylKit_gr
  
  system("nproc",intern = T) %>% as.numeric() -> nproc_system
  if(nproc_system>30){
    nproc_system = 30
  }else{
    nproc_system = nproc_system - 2
  }
  
  lapply(names(DNAm_analysis_reference_data_list$TFBS_reference)[1:4],function(nameid){
    parallel::mclapply(names(DNAm_analysis_reference_data_list$TFBS_reference[[nameid]]),mc.cores = nproc_system,function(x){
      tryCatch({
        z=DNAm_analysis_reference_data_list$TFBS_reference[[nameid]][[x]]
        y_overlap <- findOverlaps(y_gr,z)
        y_comb <- y[unique(y_overlap@from),]
        sum(y_comb$coverage*y_comb$freqC)/sum(y_comb$coverage) -> mean_beta
        sum(y_comb$coverage) -> sum_coverage
        data.frame(TF=x,beta=mean_beta,coverage=sum_coverage,regionSet=nameid)
      },error=function(e){cat ('')})
    }) %>% bind_rows() 
  }) %>% bind_rows() -> y_tfbs_result 
  return(y_tfbs_result)
}

Calculate_DNAm_for_G4AMP <- function(CpG_methylKit_filename,output_filename){
  
  require(nnls)
  require(dplyr)
  require(tidyr)
  require(parallel)
  require(readr)
  require(Matrix)
  require(matrixStats)
  require(GenomicRanges)
  
  system_name = system("hostname",intern = T)
  system_time = date()
  cat ('Running on',system_name,'at',system_time,'\n')
  x=CpG_methylKit_filename
  origin_id <- gsub('.mkdup.+','',basename(x))
  cat ('Working on',CpG_methylKit_filename,'\n')
  cat ('Reading input \n')
  readr::read_tsv(x) -> y
  y$library_name <- origin_id
  y$chr <- paste0('chr',y$chr)
  y$start = y$base %>% as.numeric()
  y$end <- y$start + 1
  y$from <- 1:nrow(y)
  y_gr <- makeGRangesFromDataFrame(y %>% dplyr::select(-base,-strand),keep.extra.columns = T)
  
  # calculate
  cat ('Working on Clock DNAm\n')
  clock_result <- Calculate_Clock_DNAm_for_G4AMP(input_methylKit_df = y,input_methylKit_gr = y_gr)
  cat ('Working on TFBS DNAm\n')
  tfbs_result <- Calculate_TFBS_DNAm_for_G4AMP(input_methylKit_df = y,input_methylKit_gr = y_gr)
  
  cat ('Saving results... \n')
  returnJsonResult <- list(
    origin = origin_id,
    file = CpG_methylKit_filename,
    clock_result = clock_result,
    tfbs_result = tfbs_result
  )
  
  saveRDS(returnJsonResult,file=paste0(gsub('.rds$','',output_filename),'.rds'))
  
  system_time = date()
  cat (' Finished on',system_name,'at',system_time,'\n')
  
}

# run 

filelist <- list.files(path='/gpfs/output/20240101_MM_YF/G4_Clinical_Test_noUMI',pattern='mkdup_CpG.methylKit$',full.names = T,recursive = T)
outputnamelist <- gsub('.mkdup_CpG.methylKit','.DNAm.simple.output.rds',basename(filelist))

mclapply(1:length(filelist),mc.cores=10,function(i){
  Calculate_DNAm_for_G4AMP(filelist[i],outputnamelist[i])
})

