# ---
# Title: iCTC Project
# Description: RCA Batch Processing
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/reboot/paper_work/Data/')
#setwd('~/Data/')
set.seed(10)

# Library Function
library(dplyr)
library(tidyr)
library(stringi)
library(stringr)
#---
# Utility Functions
#---

#---
# Method:  normalization(data,method="log")
# Description: Normalization of  the dataset. 
# Return: Normalizated data of dim(genes,samples)
#---
normalization <- function(data,method="log"){
  if (method == "log")
  {
    log_counts <- log(data + 1) / log(2)
    print(paste('Normalized Data:',dim(log_counts),sep = " "))
    return(log_counts)    
  }
  else if (method == "median")
  {
    normalized_matrix = as.matrix(data)
    total_count=colSums(data)
    med_total=median(total_count)
    for(i in 1:ncol(normalized_matrix))
    {
      normalized_matrix[,i] = data[,i] * (med_total/total_count[i])
    }
    print(paste('Normalized Data:',dim(normalized_matrix),sep = " "))
    return(normalized_matrix)
  }
  else
  {
    return("error")
  }
}

# Correlation function which return the correlation matrix #
get_corr <- function(data1,data2){
  data1[data1 == 0]<-1
  data2[data2 == 0]<-1
  correlation = cor(data1,data2,method = "pearson")
  return (correlation)
}

# Z-score compuation which return the z-score matrix #  
get_Zscore <- function(data){
  mean=colMeans(data)
  std=apply(data,2,sd)
  #return (t((t(data1)-mean)/std))
  return(scale(data,scale=std,center=mean))
}
# Stouffer score compuation which return the score matrix #    
get_stouffer <- function(data1,axis){
  return(apply(data1,axis,function(x){sum(x)/sqrt(length(x))}))
}
#---
# Utility Functions
#---

# Load the dataset
# RCA Method
load('../Data/TCGA_normalize_average_814.RData')
dim(TCGA_normalize_average_814)
load('../Data/GTEx_normalize_average_814.RData')
dim(GTEx_normalize_average_814)
filter_normalized_data<-readRDS('../Data/orginal_normalized_ctc_blood_data.rds')
dim(filter_normalized_data)
# RCA Method Batch Effect
rca_data<- na.omit(filter_normalized_data[match(toupper(rownames(TCGA_normalize_average_814)),toupper(rownames(filter_normalized_data))),])
dim(rca_data)
tcg_cor <- get_corr(rca_data,TCGA_normalize_average_814)
dim(na.omit(tcg_cor^4*sign(tcg_cor)))
tcg_zscore = get_Zscore(na.omit(tcg_cor^4)*sign(tcg_cor))
dim(tcg_zscore)
gtex_cor <- get_corr(rca_data,GTEx_normalize_average_814)
dim(na.omit(gtex_cor^4)*sign(gtex_cor))
gtex_zscore = get_Zscore(na.omit(gtex_cor^4*sign(gtex_cor)))
dim(gtex_zscore)
rca_all_data<-cbind(tcg_zscore,gtex_zscore)
dim(rca_all_data)