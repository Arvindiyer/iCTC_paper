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

#library
library(Rtsne)

# Orginal Data
filter_normalized_data<-readRDS('../Data/orginal_normalized_ctc_blood_data.rds')
tsne_data <- Rtsne(t(filter_normalized_data),
                   check_duplicates = FALSE,
                   dims = 2,
                   perplexity=50,
                   verbose=TRUE,
                   max_iter = 1000,
                   pca = TRUE,
                   pca_scale=FALSE,
                   pca_center=TRUE,
                   theta=.4,
                   eta=20,
                   exaggeration_factor=2)

# MNN Data
filter_normalized_data<-readRDS('../Data/mnn_normalized_ctc_blood_data.rds')
tsne_data <- Rtsne(t(filter_normalized_data),
                   check_duplicates = FALSE,
                   dims = 2,
                   perplexity=50,
                   verbose=TRUE,
                   max_iter = 1000,
                   pca = TRUE,
                   pca_scale=FALSE,
                   pca_center=TRUE,
                   theta=.4,
                   eta=20,
                   exaggeration_factor=2)

# Seurat Data
filter_normalized_data<-readRDS('../Data/seurat_normalized_ctc_blood_data.rds')
tsne_data <- Rtsne(t(filter_normalized_data),
                   check_duplicates = FALSE,
                   dims = 2,
                   perplexity=50,
                   verbose=TRUE,
                   max_iter = 1000,
                   pca = TRUE,
                   pca_scale=FALSE,
                   pca_center=TRUE,
                   theta=.4,
                   eta=20,
                   exaggeration_factor=2)
# RCA Data
filter_normalized_data<-readRDS('../Data/rca_normalized_ctc_blood_data.rds')
tsne_data <- Rtsne(t(filter_normalized_data),
                   check_duplicates = FALSE,
                   dims = 2,
                   perplexity=50,
                   verbose=TRUE,
                   max_iter = 1000,
                   pca = TRUE,
                   pca_scale=FALSE,
                   pca_center=TRUE,
                   theta=.4,
                   eta=20,
                   exaggeration_factor=2)
