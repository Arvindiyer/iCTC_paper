# ---
# Title: iCTC Project
# Description: Data processing Script
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# version: 0.2
# Updates:
# 1] Removed overlapping study data
# 2] Chnaged the cutoff and no-need to resample the blood data
# 3] Solved few bugs and refactor the code
# ---

#setwed 
setwd('~/final_github_iCTC/')

#load utilty functions
source('scripts/functions.R')

# Read the dataset
data <- readRDS('dataset/final_all_ctc_blood_raw_count_v2.rds')

# Dimension ofthe data
dim(data)

# Grep the total CTC and Blood
length(grep("Ctc_", colnames(data)))
length(grep("Blood_", colnames(data)))

# Apply cell filter
cell_filter_data <- cell_filter(data,10)
dim(cell_filter_data)

# Check the count of CTC and Blood
length(grep("Ctc_", colnames(cell_filter_data)))
length(grep("Blood_", colnames(cell_filter_data)))

# Apply gene filtering to the dataset
gene_filter_data <- gene_filter(cell_filter_data,5,10)
dim(gene_filter_data)

CTC_GSE51827_sample=grep("Ctc_GSE51827", colnames(gene_filter_data))
CTC_GSE55807_sample=grep("Ctc_GSE55807", colnames(gene_filter_data))
CTC_GSE60407_sample=grep("Ctc_GSE60407", colnames(gene_filter_data))
CTC_GSE67939_sample=grep("Ctc_GSE67939", colnames(gene_filter_data))
CTC_GSE67980_sample=grep("Ctc_GSE67980", colnames(gene_filter_data))
CTC_GSE74639_sample=grep("Ctc_GSE74639", colnames(gene_filter_data))
CTC_GSE75367_sample=grep("Ctc_GSE75367", colnames(gene_filter_data))
CTC_GSE109761_sample = grep("Ctc_GSE109761",colnames(gene_filter_data))
CTC_GSE86978_sample = grep("Ctc_GSE86978",colnames(gene_filter_data))
CTC_GSE38495_sample = grep("Ctc_GSE38495",colnames(gene_filter_data))
CTC_Naveen_sample = grep('Ctc_Naveen',colnames(gene_filter_data))
Blood_Satijalab_sample=grep("Blood_Satijalab", colnames(gene_filter_data))
Blood_PBMC_5419_sample=grep("Blood_PBMC_5419", colnames(gene_filter_data))
Blood_GSE67980_sample=grep("Blood_GSE67980", colnames(gene_filter_data))
Blood_GSE109761_sample=grep("Blood_GSE109761", colnames(gene_filter_data))
Blood_GSE81861_sample=grep("Blood_GSE81861", colnames(gene_filter_data))
Blood_Genomics_sample=grep("Blood_Genomics", colnames(gene_filter_data))

index <- c(CTC_GSE51827_sample,
           CTC_GSE55807_sample,
           CTC_GSE60407_sample,
           CTC_GSE67939_sample,
           CTC_GSE67980_sample,
           CTC_GSE74639_sample,
           CTC_GSE75367_sample,
           CTC_GSE109761_sample,
           CTC_GSE86978_sample,
           CTC_GSE38495_sample,
           CTC_Naveen_sample,
           Blood_Satijalab_sample,
           Blood_PBMC_5419_sample,
           Blood_GSE67980_sample,
           Blood_GSE109761_sample,
           Blood_GSE81861_sample,
           Blood_Genomics_sample)

# Total samples 
length(index)
gene_filter_data_reduced <- gene_filter_data[,index]
dim(gene_filter_data_reduced)

# Median Normalization
normalization_data <- normalization(gene_filter_data_reduced,method = "median")
saveRDS(normalization_data,file='dataset/final_all_ctc_blood_median_normalized_data_v2.rds')

# Machine Learning Gene Selection
# Read marker genes
genes <- read_csv('dataset/experimental_marker_genes.csv')
ml_data <- normalization_data[intersect(rownames(normalization_data),toupper(genes$genes)),]
dim(ml_data)
saveRDS(ml_data,file='dataset/ml_median_normalized_data_v2.rds')
ml_genes <- intersect(rownames(normalization_data),toupper(genes$genes))