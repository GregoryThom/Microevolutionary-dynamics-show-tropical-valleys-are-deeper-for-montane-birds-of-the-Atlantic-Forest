library(PipeMaster)
library(ggplot2)
#library(caret)
#library(doMC)
library(ggpubr)
setwd("../inputs_and_data_sets/7_loci_fasta/loci_all_sp/1_pipe_master_10K/")

#stat_aurulentus <- read.table("../dados_espaciais/summary_stats/aurulentus/1_aurulentus_sum_stat", header = T)


###drawing a model so Pipemaster can calculate the summary stats
#m1 isolation
m1 <- main.menu()
###m8 isolation migration

files <- gsub(".pops", "", list.files(path="../1_pipe_master_10K/", pattern=".pops", full.names=F, recursive=FALSE))
file <- files[[1]]
for(file in files){
  setwd("../inputs_and_data_sets/7_loci_fasta/loci_all_sp/1_pipe_master_10K/")
  
  pop.assign <- read.delim(paste(file, ".pops", sep = ""), header = FALSE)
  pop.assign$V2 <- 1 
  
  m1 <- get.data.structure(m1,path.to.fasta = paste("../inputs_and_data_sets/7_loci_fasta/loci_all_sp/1_pipe_master_10K/", file, sep = ""), pop.assign=pop.assign)
  
  
  ##observed_sum_stats
  
  stat <- obs.sumstat.ngs(model = m1, path.to.fasta = paste("../inputs_and_data_sets/7_loci_fasta/loci_all_sp/1_pipe_master_10K/", file, sep = ""), 
                          pop.assign = pop.assign, moments = F)
  
  setwd("../inputs_and_data_sets/7_loci_fasta/loci_all_sp/1_pipe_master_10K/")
  write.table(stat, file = paste("1_", file, "_sum_stat",sep = ""))
  
}