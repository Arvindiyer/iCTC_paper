# ---
# Title: iCTC Project
# Description: Main Script to do the Network Simulation
# Contact: Kishore Hari <kishorehari@iisc.ac.in>
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in>, Kishore Hari <kishorehari@iisc.ac.in>, Mohit Kumar Jolly <mkjolly@iisc.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

library(tidyverse)
library(magrittr)
library(stringr)
library(ggthemes)
library(parallel)
wd <-  getwd()
racipe <- paste0(wd, "/RACIPE-1.0-master/")

Source <- c("ZEB", "miR200", "ZEB", "SNAIL", "SNAIL", "GRHL2", "ZEB", "OVOL2", "ZEB", "GRHL2", 
            "ZEB", "STEP1", "CLDN4", "OVOL2", "GRHL2", "ZEB", "STEP1", "miR200", "STAT3")
Target <- c("miR200", "ZEB", "ZEB", "miR200", "ZEB", "ZEB", "GRHL2", "ZEB", "OVOL2", "OVOL2", 
            "CDH1", "CDH1", "CDH1", "CDH1", "CDH1", "VIM1", "VIM1", "STAT3", "ZEB")
Type <- c(2, 2, 1, 2, 1, 2,2,2,2,1, 2, 1, 1, 1, 1, 1,2, 2, 1)

df <- data.frame(Source = Source, Target = Target, Type = Type)

write_delim(df, "network.topo", delim = "\t")

setwd(wd)
for (i in 1:3)
{
    if (!dir.exists(as.character(i))) dir.create(as.character(i))
}

system(paste0("cp network.topo ", racipe))
setwd(racipe)

cl <- makeCluster(3)
clusterExport(cl, "wd")

parLapplyLB(cl, 1:3, function(i){
    system(paste0("cp network.topo network_",i,".topo"))
    system(paste0("./RACIPE network_",i,".topo -num_paras 10000 > network_",i,".txt"))
    system(paste0("cp network_",i,"* ", wd, i))
})

stopCluster(cl)

setwd(wd)
source(paste0(wd, "/ss_and_phase_analysis.R"))



    sapply(as.character(1:3), function(y){
        f <- list.files(path = y, pattern = "nfull")
        #if(length(f) == 0)
            z_score_normalized_calculator(y)
    })


    sapply(as.character(1:3), function(y){
        f <- list.files(path = y, pattern = "lz")
        #if(length(f) == 0)
            z_score_calculator_log(y)
    })

  setwd("..")
  steady_state_patterns(folder = "./IIITD", n_cluster = 2, scatter = c(4,8))
  steady_state_patterns(folder = "./IIITD", n_cluster = 3, scatter = c(4,8))
  frequency_plot_maker(folder = "./IIITD", n_cluster = 2, primary_comps = c(1,2))
  frequency_plot_maker(folder = "./IIITD", n_cluster = 3, primary_comps = c(4,8))
  Quadrant_plots(folders = "./IIITD", primary = c(1,2))
  #Phase_plots(folders = "./IIITD", primary = c(1,2), short_code = T)
  #boolean_converter("./IIITD", 10000)
  
