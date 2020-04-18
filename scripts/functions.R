# ---
# Title: iCTC Project
# Description: Utility Functions
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# version: 0.2
# Updates:
# 1] Added hyrdo-seq data
# 2] Solved few bugs and refactor the code
# ---

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
    log_counts <- log2(data + 1)
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
# Method:  dispersion_genes(data,ngenes_keep)
# Description: Feature selection method.
# Return: Feature Selected Genes based on coefficent of variation 
# Credits: Code snippet taken from DropClust Package
#---
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
#---
# Method:  get_Zscore(data)
# Description:  Z-score compuation method
# Return: Z-score matrix 
#---
get_Zscore <- function(data){
  data[data == 0]<-1
  mean=colMeans(data)
  std=apply(data,2,sd)
  #return (t((t(data1)-mean)/std))
  return(scale(data,scale=std,center=mean))
}
#---
# Method:  get_stouffer(data,axis)
# Description:  Return stouffer or phenotype score
# Params: data: z-score matrix axis: computation to be done on row or column
# Return: Z-score matrix 
#---
get_stouffer <- function(data1,axis){
  return(apply(data1,axis,function(x){sum(x)/sqrt(length(x))}))
}
#---
# Method:  round_df(data,digits)
# Description:  Round of dataframe
# Return: Rounded dataframe
# Credits: Stackoverflow
#---
round_df <- function(df, digits = 3) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}
#---
# Method:  round_df(data,digits)
# Description:  Moving average smoother for columns of a matrix
# Return: Rounded dataframe
# Credits: Taken from edgeR Package by Gordon Smyth
#---
movingAverageByCol <- function(x,width=5,full.length=TRUE)
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
# Median for box plot
fun_mean <- function(x){
  return(data.frame(y=median(x),label=median(x,na.rm=T)))
}
# Function to conver the count data to tpm
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
# ---
# Utility Functions End
# ---