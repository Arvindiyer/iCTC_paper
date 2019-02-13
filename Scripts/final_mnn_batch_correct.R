# ---
# Title: iCTC Project
# Description: MNNCorrect Batch Processing
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
library(scran)

# Load the dataset
normalized_matrix_filter<-readRDS('../Data/orginal_normalized_ctc_blood_data.rds')
dim(normalized_matrix_filter)
# MNN Method Batch Effect
CTC_GSE51827_sample=grep("CTC_GSE51827", colnames(normalized_matrix_filter))
CTC_GSE55807_sample=grep("CTC_GSE55807", colnames(normalized_matrix_filter))
CTC_GSE60407_sample=grep("CTC_GSE60407", colnames(normalized_matrix_filter))
CTC_GSE67939_sample=grep("CTC_GSE67939", colnames(normalized_matrix_filter))
CTC_GSE67980_sample=grep("CTC_GSE67980", colnames(normalized_matrix_filter))
CTC_GSE74639_sample=grep("CTC_GSE74639", colnames(normalized_matrix_filter))
CTC_GSE75367_sample=grep("CTC_GSE75367", colnames(normalized_matrix_filter))
Blood_GSE81861_sample=grep("Blood_GSE81861", colnames(normalized_matrix_filter))
Blood_Satijalab_sample=grep("Blood_Satijalab", colnames(normalized_matrix_filter))
Blood_PBMC_5419_sample=grep("Blood_PBMC_5419", colnames(normalized_matrix_filter))
Blood_genomics_10x_mtx=grep("Blood_genomics_10x_mtx", colnames(normalized_matrix_filter))
Blood_GSE67980_sample=grep("Blood_GSE67980", colnames(normalized_matrix_filter))

# All Data
mnn_data<-mnnCorrect(
  log(1+normalized_matrix_filter[,CTC_GSE51827_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE55807_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE60407_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE67939_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE67980_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE74639_sample]),  
  log(1+normalized_matrix_filter[,CTC_GSE75367_sample]),  
  log(1+normalized_matrix_filter[,Blood_GSE81861_sample]),
  log(1+normalized_matrix_filter[,Blood_Satijalab_sample]),
  log(1+normalized_matrix_filter[,Blood_PBMC_5419_sample]),
  log(1+normalized_matrix_filter[,Blood_genomics_10x_mtx]),
  log(1+normalized_matrix_filter[,Blood_GSE67980_sample]),
  k=5
)
final_run_mnn <- cbind(mnn_data$corrected[[1]],
                       mnn_data$corrected[[2]],
                       mnn_data$corrected[[3]],
                       mnn_data$corrected[[4]],
                       mnn_data$corrected[[5]],
                       mnn_data$corrected[[6]],
                       mnn_data$corrected[[7]],
                       mnn_data$corrected[[8]],
                       mnn_data$corrected[[9]],
                       mnn_data$corrected[[10]],
                       mnn_data$corrected[[11]],
                       mnn_data$corrected[[12]])
dim(final_run_mnn)

colnames(final_run_mnn) <- colnames(cbind(  normalized_matrix_filter[,CTC_GSE51827_sample],
                                            normalized_matrix_filter[,CTC_GSE55807_sample],
                                            normalized_matrix_filter[,CTC_GSE60407_sample],
                                            normalized_matrix_filter[,CTC_GSE67939_sample],
                                            normalized_matrix_filter[,CTC_GSE67980_sample],
                                            normalized_matrix_filter[,CTC_GSE74639_sample],  
                                            normalized_matrix_filter[,CTC_GSE75367_sample],  
                                            normalized_matrix_filter[,Blood_GSE81861_sample],
                                            normalized_matrix_filter[,Blood_Satijalab_sample],
                                            normalized_matrix_filter[,Blood_PBMC_5419_sample],
                                            normalized_matrix_filter[,Blood_genomics_10x_mtx],
                                            normalized_matrix_filter[,Blood_GSE67980_sample]))
# Only CTC
mnn_data<-mnnCorrect(
  log(1+normalized_matrix_filter[,CTC_GSE51827_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE55807_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE60407_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE67939_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE67980_sample]),
  log(1+normalized_matrix_filter[,CTC_GSE74639_sample]),  
  log(1+normalized_matrix_filter[,CTC_GSE75367_sample]),
  k=5
)
final_run_mnn <- cbind(mnn_data$corrected[[1]],
                       mnn_data$corrected[[2]],
                       mnn_data$corrected[[3]],
                       mnn_data$corrected[[4]],
                       mnn_data$corrected[[5]],
                       mnn_data$corrected[[6]],
                       mnn_data$corrected[[7]])

colnames(final_run_mnn) <- colnames(cbind(  normalized_matrix_filter[,CTC_GSE51827_sample],
                                            normalized_matrix_filter[,CTC_GSE55807_sample],
                                            normalized_matrix_filter[,CTC_GSE60407_sample],
                                            normalized_matrix_filter[,CTC_GSE67939_sample],
                                            normalized_matrix_filter[,CTC_GSE67980_sample],
                                            normalized_matrix_filter[,CTC_GSE74639_sample],  
                                            normalized_matrix_filter[,CTC_GSE75367_sample]))  