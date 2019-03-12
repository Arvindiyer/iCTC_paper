# ---
# Title: iCTC Project
# Description: Survival analysis of E_M combination on TCGA Data
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Load utility functions
source('scripts/functions.R')
# Library
library(survival)
library(survminer)
library(survcomp)
library(randomForestSRC)

#---
# GBM
#---

# Code used to do KM analyis.

# # Load the epi and mesen genes
# epi_genes <- readRDS('dataset/new_final_filter_epitheilal_genes.rds')
# mesn_genes <- readRDS('dataset/new_final_filter_mesenchymal_genes.rds')
# expression_data <- readRDS('dataset/final_marker_normalized_matrix.rds')
# load('dataset/gbm_final_mRNA_raw_count_epi_mesn.data')
# duplicated.columns <- duplicated(t(gbm_mRNA))
# gbm_mRNA <- gbm_mRNA[, !duplicated.columns]
# dim(gbm_mRNA)
# mRNA <- normalization(gbm_mRNA,method = "median")
# #mRNA<-t(mRNA)
# dim(mRNA)
# 
# # Load Clincial Data
# load('dataset/gbm_final_clincial.data')
# #dim(unique(gbm_clincial))
# gbm_clincial<-unique(gbm_clincial)
# #gbm_clincial$gender <- factor(gbm_clincial$gender)
# gbm_clincial[gbm_clincial$gender=='male',]$gender<-0
# gbm_clincial[gbm_clincial$gender=='female',]$gender<-1
# gbm_clincial$gender <- as.numeric(as.character(gbm_clincial$gender))
# str(gbm_clincial)
# # Load Survival
# load('dataset/gbm_final_survival.data')
# #dim(unique(gbm_clincial))
# gbm_survival<-unique(gbm_survival)
# #gbm_survival$status<-factor(gbm_survival$status)
# gbm_survival[gbm_survival$status=='alive',]$status<-0
# gbm_survival[gbm_survival$status=='dead',]$status<-1
# #gbm_survival$time<-gbm_survival$time
# #gbm_survival[gbm_survival$time<0,]$time=0
# gbm_survival$status <- as.numeric(as.character(gbm_survival$status))
# str(gbm_survival)
# # Select one of the combination  
# comb <- combn(epi_genes,2,simplify = F)
# comb <- combn(mesn_genes,2,simplify = F)
# comb <- expand.grid(epi_genes, mesn_genes,stringsAsFactors = FALSE)
# str(comb)
# length(comb)
# dim(comb)
# p_val <- list()
# #z_score_mRNA <- get_Zscore(log(mRNA)+1)
# #dim(z_score_mRNA)
# # Run for loop dependng on comb object and comment appropiate t inside the for loop
# for (val in 1:nrow(comb)) {
#   #for (val in comb) {
#   t<-c(comb[val,]$Var1,comb[val,]$Var2)
#   #  t<-c(val[1],val[2])
#   gbm_exp <- mRNA[t,]
#   gbm_zscore <- get_Zscore(t(log(gbm_exp+1)))
#   gbm_stouffer <- get_stouffer(gbm_zscore,1)
#   med<-median(gbm_stouffer)
#   cluster=sapply(gbm_stouffer,function(x){if (x >= med) return("High") else return("Low")})
#   surv_data <- data.frame('time'=gbm_survival$time,'event'=gbm_survival$status,'cluster'=cluster)
#   km <- surv_fit(Surv(time, event) ~ cluster,data = surv_data)
#   p_val[paste(t[1],t[2],sep = "_")]<-surv_pvalue(km)$pval
#   # plot<-ggsurvplot(km,
#   #                  pval = TRUE, 
#   #                  title = paste("KM Plot for combination",t[1],t[2]),
#   #                  risk.table = TRUE, # Add risk table
#   #                  risk.table.col = "strata", # Change risk table color by groups
#   #                  linetype = "strata", # Change line type by groups
#   #                  surv.median.line = "hv", # Specify median survival
#   #                  #ggtheme = theme_bw(), # Change ggplot2 theme
#   #                  palette = "Set1")
#   # ggsave(filename = paste("epi_mesn",t[1],t[2],".pdf",sep = "_"),
#   #        plot = print(plot,newpage = FALSE),
#   #        path = "../../new_survival/")
#   print(paste("new_epi_mesn",t[1],t[2],"Done",sep = "_"))
# }

# epi_combination_p_val <- unlist(p_val)
# length(epi_combination_p_val)
# save(epi_combination_p_val,file = "dataset/gbm/final_updated_new_epi_combination_pval.Rdata")
# 
# mesn_combination_p_val <- unlist(p_val)
# length(mesn_combination_p_val)
# save(mesn_combination_p_val,file = "dataset/gbm/final_updated_new_mesn_combination_pval.Rdata")
# 
# epi_mesn_combination_p_val <- unlist(p_val)
# length(epi_mesn_combination_p_val)
# save(epi_mesn_combination_p_val,file = "dataset/gbm/final_updated_new_epi_mesn_combination_pval.Rdata")


load('dataset/gbm/final_updated_new_epi_combination_pval.Rdata')
load('dataset/gbm/final_updated_new_mesn_combination_pval.Rdata')
load('dataset/gbm/final_updated_new_epi_mesn_combination_pval.Rdata')

p_val_data <- c(epi_combination_p_val,mesn_combination_p_val,epi_mesn_combination_p_val)
label<-c(rep('Epithelial',120),rep('Mesenchymal',903),rep('Combination',688))
data<-data.frame('p values'= -log10(p_val_data),'method'=label)
dim(data)

ggboxplot(data,
          x = "method",
          y = "p.values",
          color = "method",
          bxp.errorbar = TRUE,
          notch = FALSE,
          palette = "Set2",
          outlier.shape = NA,
          add = c("jitter"),
          add.params =list("size"= .75)
)+
  guides(size = FALSE)+
  stat_summary(fun.y = median, geom="point",colour="darkred", size=3)+
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7,aes( label=round(..y.., digits=3)))+
  theme(axis.text.y = element_text(family = "Helvitca",colour = "black",size = 12),
        axis.text.x = element_text(family = "Helvitca",colour = "black",size = 12,angle = 0, hjust = 0.5),
        legend.position = "None",
        legend.text = element_text(family = "Helvitca",color = "black",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"))+
  xlab('')+
  #coord_cartesian(ylim = c(0,1.75))+
  ylab('-Log10 p-values')

combine_data<-data.frame('p values'= p_val_data,'method'=label)
reduced_data <- combine_data[combine_data$p.values < 0.05,]

table(combine_data$method)
table(reduced_data$method)
View(reduced_data)
table(reduced_data$method)*100/table(combine_data$method)
View(reduced_data)

# Similar way for other 3 cancer