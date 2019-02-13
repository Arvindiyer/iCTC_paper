# ---
# Title: iCTC Project
# Description: Supplementary Figure-6
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/reboot/paper_work/Data/')
#setwd('~/Data/')
set.seed(10)

library(ggpubr)
library(ggplot2)
library(RColorBrewer)

# ---
# Utility Functions
# ---

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

dispersion_genes<-function(normalized_data, ngenes_keep=1000) {
  avg<-Matrix::colMeans(normalized_data$m)
  vars <- apply(normalized_data$m,2,var)
  dispersion = vars/avg
  # ColDispersion(normalized_data$m)
  avg_bin<-cut(avg,breaks=c(-Inf,quantile(avg,seq(0.1,1,0.05)),Inf))
  df=as.data.frame(cbind(avg,dispersion,avg_bin))
  
  variance_by_bin<-plyr::ddply(df,"avg_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-variance_by_bin$bin_median[match(df$avg_bin,variance_by_bin$avg_bin)]
  df$bin_disp_mad<-variance_by_bin$bin_mad[match(df$avg_bin,variance_by_bin$avg_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  
  print("Sort Top Genes...")
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]
  print("Cutoff Genes...")
  df$used<-df$dispersion_norm >= disp_cut_off
  features = head(order(-df$dispersion_norm),ngenes_keep)
  dispersed_genes = normalized_data$use_genes[features]
  
  return(dispersed_genes)
}
# Z-score compuation which return the z-score matrix #  
get_Zscore <- function(data){
  data[data == 0]<-1
  mean=colMeans(data)
  std=apply(data,2,sd)
  #return (t((t(data1)-mean)/std))
  return(scale(data,scale=std,center=mean))
}
# Stouffer score compuation which return the score matrix #    
get_stouffer <- function(data1,axis){
  return(apply(data1,axis,function(x){sum(x)/sqrt(length(x))}))
}

round_df <- function(df, digits = 3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

movingAverageByCol <- function(x,width=5,full.length=TRUE)
  #	Moving average smoother for columns of a matrix
  #	Gordon Smyth
  #	17 Feb 2011
  # Taken from edgeR package
{
  x <- as.matrix(x)
  width <- as.integer(width)
  if(width<=1) return(x)
  n <- nrow(x)
  m <- ncol(x)
  if(width>n) {
    width <- n
    warning("reducing moving average width to nrow(x)")
  }
  if(full.length) {
    half1 <- ceiling(width/2)
    half2 <- floor(width/2)
    x <- rbind(matrix(0,half1,m),x,matrix(0,half2,m))
  } else {
    if(width==n) return(matrix(colMeans(x),1L,m))
    x <- rbind(matrix(0,1,m),x)
  }
  n2 <- nrow(x)
  x <- apply(x,2,cumsum)
  x <- x[(width+1):n2,,drop=FALSE]-x[1:(n2-width),,drop=FALSE]
  n3 <- nrow(x)
  w <- rep(width,n3)
  if(full.length) {
    if(half1>1) w[1:(half1-1)] <- width-(half1-1):1
    w[(n3-half2+1):n3] <- width-(1:half2)
  }
  x/w
}

# ---
# Utility Functions End
# ---

# data <- readRDS('ctc_data.rds')
# raw.counts<-apply(data,2,function(x) {storage.mode(x) <- 'integer'; x})
# rm(data)
# dim(raw.counts)
# # Apply cell filtering to the dataset
# cell_filter_data <- cell_filter(raw.counts,percentage = 5)
# dim(cell_filter_data)
# # Apply gene filtering to the dataset
# gene_filter_data <- gene_filter(cell_filter_data,5,92)
# dim(gene_filter_data)
# # Apply Median Normalization
# normalized_data<-normalization(gene_filter_data,"median")
# dim(normalized_data)

# Loading Genes
filter_epi_genes <- readRDS('final_filter_epitheilal_genes.rds')
filter_mesn_genes <- readRDS('final_filter_mesenchymal_genes.rds')
filter_csc_genes <- readRDS('final_filter_cancer_stem_genes.rds')
normalized_data <- readRDS('final_marker_normalized_matrix.rds')
new_mat_to_plot_gene<-unique(normalized_data)
dim(new_mat_to_plot_gene)
length(na.omit(match(unique(filter_epi_genes),rownames(new_mat_to_plot_gene))))
length(na.omit(match(unique(filter_mesn_genes),rownames(new_mat_to_plot_gene))))
length(na.omit(match(unique(filter_csc_genes),rownames(new_mat_to_plot_gene))))

CTC_GSE51827_sample=length(grep("CTC_GSE51827", colnames(new_mat_to_plot_gene)))
CTC_GSE55807_sample=length(grep("CTC_GSE55807", colnames(new_mat_to_plot_gene)))
CTC_GSE60407_sample=length(grep("CTC_GSE60407", colnames(new_mat_to_plot_gene)))
CTC_GSE67939_sample=length(grep("CTC_GSE67939", colnames(new_mat_to_plot_gene)))
CTC_GSE67980_sample=length(grep("CTC_GSE67980", colnames(new_mat_to_plot_gene)))
CTC_GSE74639_sample=length(grep("CTC_GSE74639", colnames(new_mat_to_plot_gene)))
CTC_GSE75367_sample=length(grep("CTC_GSE75367", colnames(new_mat_to_plot_gene)))

col_label<-c(rep("CTC Breast #1 (Aceto N et al.)",CTC_GSE51827_sample),
             rep("CTC Breast #2 (Yu M et al.)",CTC_GSE55807_sample),
             rep("CTC Pancreas #3 (Ting DT et al.)",CTC_GSE60407_sample),
             rep("CTC Breast #4 (Sarioglu AF et al.)",CTC_GSE67939_sample),
             rep("CTC Prostrate #5 (Miyamoto DT et al.)",CTC_GSE67980_sample),
             rep("CTC Lung #6 (Zheng Y et al.)",CTC_GSE74639_sample),
             rep("CTC Breast #7 (Jordan NV et al.)",CTC_GSE75367_sample))

row_annotation = c(rep('Epithelial',11),rep('Mesenchymal',25),rep('Cancer Stem Cell',15))
length(row_annotation)
length(col_label)
dim(new_mat_to_plot_gene)
z_score <- get_Zscore(log(new_mat_to_plot_gene+1))
dim(na.omit(z_score))
dim(unique(z_score))
annotation_col <- data.frame(Study = factor(col_label))
rownames(annotation_col) <- colnames(z_score)
dim(annotation_col)
annotation_rows <- data.frame(Marker = factor(row_annotation))
rownames(annotation_rows)<-rownames(new_mat_to_plot_gene)
dim(annotation_rows)
test1<-get_stouffer(z_score[1:11,],2)
test2<-get_stouffer(z_score[12:36,],2)
test3<-get_stouffer(z_score[37:51,],2)
new_data <- data.frame('Cancer'=test3,'Epitheial'=test1,'Mesenchymal'=test2)
dim(new_data)
ann_colors = list(Study = c("CTC Breast #1 (Aceto N et al.)"="#E41A1C",
                            "CTC Breast #2 (Yu M et al.)"="#377EB8",
                            "CTC Pancreas #3 (Ting DT et al.)"="#4DAF4A",
                            "CTC Breast #4 (Sarioglu AF et al.)"="#984EA3",
                            "CTC Prostrate #5 (Miyamoto DT et al.)"="#FF7F00",
                            "CTC Lung #6 (Zheng Y et al.)"="#FFFF33",
                            "CTC Breast #7 (Jordan NV et al.)"="#A65628"))

new_data$Study <- col_label
dim(new_data)

matrix_data <- new_data
E<- matrix_data$Epitheial+abs(min(matrix_data$Epitheial))+1
M<- matrix_data$Mesenchymal+abs(min(matrix_data$Mesenchymal))+1
E_M <- E/M
matrix_data$E_M <-E_M
dim(matrix_data)
order_matrix_data <- matrix_data[order(matrix_data$E_M),]
dim(order_matrix_data)
head(order_matrix_data)
# Scatter Plot
ggplot(round_df(order_matrix_data,2), aes(x=Epitheial, y=Mesenchymal,color=Study)) +
  facet_wrap(. ~ Study) +
  scale_color_brewer(palette="Set1")+
  geom_point() +
  theme_bw()+
  theme(axis.text = element_text(family = "Helvitca",size = 10),
        legend.position = "right",
        legend.title = element_blank(),
        strip.text = element_text(family = "Helvitca",size = 5),
        legend.text = element_text(family = "Helvitca",size = 8))+
  guides(color = guide_legend(override.aes = list(size = 2)))+
  xlab('Epithelial Singature')+
  ylab('Mesenchymal Singature')


#write.csv(order_matrix_data,file="E_M_order_data.csv")
