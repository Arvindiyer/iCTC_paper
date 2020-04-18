# ---
# Title: iCTC Project
# Description: Differntial Expression Analysis
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

library(superheat)
library(RColorBrewer)
source('scripts/functions.R')

# Load the dataset.
normalized_data <- readRDS('dataset/final_all_ctc_blood_median_normalized_data_v2.rds')
dim(normalized_data)
length(grep("Ctc_", colnames(normalized_data)))

# Perform Wilcoxon test on ctc vs blood and compute pvalues
pval<-apply(normalized_data,1,function(x) wilcox.test(x[1:538],x[538:1861])$p.value)
# Perform Bonferroni correction on p value to get q value.
fdr <- p.adjust(pval,method="fdr")
# Prepare Differnetial final data with fold change value
DE_res <- data.frame(cbind(pvalues=pval,qvalues=fdr))
FC = log2(Matrix::rowMeans(normalized_data[,1:538])/Matrix::rowMeans(normalized_data[,538:1861]))
df<-data.frame(gene = rownames(DE_res),p_val=DE_res$pvalues,q_val = DE_res$qvalues,fc = FC)
dim(df)

# Save the data
save(df,file='dataset/new_de_data.Rdata')

# Load the diffenetial data
load('dataset/new_de_data.Rdata')
reduced_df<- df[df$q_val< 0.05,]
dim(reduced_df)
genes<-read.csv('data/cell_surface_marker.csv',sep = "\t",header = TRUE)
common <- intersect(toupper(rownames(df)),toupper(genes$Gene))
new_Data <- normalized_data[common,]
# Picking only thise genes which are expressed in atleast 589 CTC samples i.e 80% of the data
express <- apply(new_Data[,1:538],1,function(x) sum(x>0)>=(430))
sum(express)
new_Data_filter_ctc_genes <- rownames(new_Data[express,])
length(new_Data_filter_ctc_genes)
data_df<-na.omit(reduced_df[new_Data_filter_ctc_genes,])
dim(data_df)
data_df<-data_df[order(data_df$fc),]
# Generate Average Expression for CTC
dim(normalized_data[rownames(data_df),])
expression_data<-normalized_data[rownames(data_df),1:538]
rowMeans(expression_data)
View(rowMeans(expression_data))
# Generate Average Expression for Blood
expression_data<-normalized_data[rownames(data_df),538:1861]
rowMeans(expression_data)
View(rowMeans(expression_data))