# ---
# Title: iCTC Project
# Description: Raw Data Processing
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/reboot/paper_work/Data/')
#setwd('~/Data/')
set.seed(10)

#---
# Utility Functions
#--

#---
# Method:  cell_filter(data,percentage)
# Description: Cell filtering in the dataset e.x Keep only those cells in which 5% of genes are expressed.
# Return: Filter data of dim(genes,samples)
#---
cell_filter <- function(data,percentage=10){
  express=apply(data,2,function(x) sum(x>0)>=(length(data[,1])*(percentage/100)))
  print(paste('Filtered Samples:',sum(express),sep = " "))
  return(data[,express])
}
#---
# Method:  gene_filter(data,min.count=2, min.cell=3)
# Description: Gene filtering in the dataset e.x Genes with count greater than min.count = 2 in atleast min.cell = 3 cells is retained.
# Return: Filter data of dim(genes,samples)
#---
gene_filter <- function(data,min.count=2,min.cell=3){
  express=apply(data,1,function(x) sum(x>min.count)>=(min.cell))
  print(paste('Filtered Genes:',sum(express),sep = " "))
  return(data[express,])
}

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
#---
# Utility Functions End
#--

# Load the data
load('../Data/ctc_blood_data.Rdata')
dim(data)
raw.counts<-apply(data,2,function(x) {storage.mode(x) <- 'integer'; x})
rm(data)
dim(raw.counts)
# Apply cell filtering to the dataset
cell_filter_data <- cell_filter(raw.counts,percentage = 6)
dim(cell_filter_data)
# Check the count of CTC and Blood
length(grep("CTC_", colnames(cell_filter_data)))
length(grep("Blood_", colnames(cell_filter_data)))
# Apply gene filtering to the dataset
gene_filter_data <- gene_filter(cell_filter_data,5,10)
dim(gene_filter_data)
# Apply median Normalization
normalized_data<- normalization(gene_filter_data,method = "median")
dim(normalized_data)
saveRDS(normalized_data,file='ctc_blood_normalized_data.rds')
# Gene Selection based on markers genes curaed from marker gene database.
genes<-read.csv('experimental_marker_genes.csv')
dim(normalized_data)
index<- na.omit(match(toupper(genes$genes),rownames(normalized_data)))
length(index)
final_data <- normalized_data[index,]
dim(final_data)
#saveRDS('../Data/orginal_normalized_ctc_blood_data.rds')