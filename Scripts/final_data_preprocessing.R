# ---
# Title: iCTC Project
# Description: Raw Data and Preprocessing
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/reboot/paper_work/Data/')
#setwd('~/Data/')
set.seed(10)


# Use readRDS to read the datasets of which  GSE67980 study consist both blood and ctc data
GSE51827<- readRDS(file="GSE51827_CTC.rds")
GSE55807<- readRDS(file="GSE55807_CTC.rds")
GSE60407<- readRDS(file="GSE60407_CTC.rds")
GSE67939<- readRDS(file="GSE67939_CTC.rds")
GSE67980<- readRDS(file="GSE67980_CTC.rds")
GSE74639<- readRDS(file="GSE74639_CTC.rds")
GSE75367<- readRDS(file="GSE75367_CTC.rds")
pbmc_3k<- readRDS(file="pbmc_3k_blood.rds")
GSE81861<- readRDS(file="GSE81861_blood.rds")
pbmc_6k<- readRDS(file="pbmc_6k_blood.rds")
EGAS00001002560<- readRDS(file="EGAS00001002560_blood.rds")


# Find the common genes among these studies.
intersect_genes <- Reduce(intersect, list(rownames(GSE51827),
                                          rownames(GSE55807),
                                          rownames(GSE60407),
                                          rownames(GSE67939),
                                          rownames(GSE67980),
                                          rownames(GSE74639),
                                          rownames(GSE75367),
                                          rownames(pbmc_3k),
                                          rownames(GSE81861),
                                          rownames(pbmc_6k),
                                          rownames(EGAS00001002560)))
length(intersect_genes)
#saveRDS(intersect_genes,file="common_ctc_blood_genes.rds")
#write.csv(intersect_genes,file="common_ctc_blood_genes.csv")
# Saving the comibined data
data <- cbind(  GSE51827[intersect_genes,],
                GSE55807[intersect_genes,],
                GSE60407[intersect_genes,],
                GSE67939[intersect_genes,1:15],
                GSE67980[intersect_genes,1:164],
                GSE74639[intersect_genes,],
                GSE75367[intersect_genes,],
                GSE81861[intersect_genes,1:123],
                pbmc_3k[intersect_genes,],
                pbmc_6k[intersect_genes,],
                GSE67980[intersect_genes,165:169],
                as.data.frame(as.matrix(EGAS00001002560[intersect_genes,]))
                )
# Setting New Colnames
new_col_names= c(
  c(paste("CTC_GSE51827", colnames(GSE51827) ,c(1:29),sep="_"))
  ,c(paste("CTC_GSE55807", colnames(GSE55807) ,c(1:6),sep="_"))
  ,c(paste("CTC_GSE60407", colnames(GSE60407) ,c(1:7),sep="_"))
  ,c(paste("CTC_GSE67939", colnames(GSE67939[,1:15]) ,c(1:15),sep="_"))
  ,c(paste("CTC_GSE67980", colnames(GSE67980[,1:164]) ,c(1:164),sep="_"))
  ,c(paste("CTC_GSE74639", colnames(GSE74639) ,c(1:16),sep="_"))
  ,c(paste("CTC_GSE75367", colnames(GSE75367) ,c(1:74),sep="_"))
  ,c(paste("Blood_GSE81861", colnames(GSE81861[,1:123]) ,c(1:123),sep="_"))
  ,c(paste("Blood_Satijalab", colnames(pbmc_3k) ,c(1:2700),sep="_"))
  ,c(paste("Blood_PBMC_5419", colnames(pbmc_6k) ,c(1:5419),sep="_"))
  ,c(paste("Blood_GSE67980", colnames(GSE67980[,165:169]) ,c(1:5),sep="_"))
  ,c(paste("Blood_genomics_10x_mtx", colnames(EGAS00001002560) ,c(1:28855),sep="_"))
)
colnames(data)<-new_col_names
dim(data)

#saveRDS(data,file='ctc_blood_data.rds')
# Find the common genes among only CTC studies
intersect_genes <- Reduce(intersect, list(rownames(GSE51827),
                                          rownames(GSE55807),
                                          rownames(GSE60407),
                                          rownames(GSE67939),
                                          rownames(GSE67980),
                                          rownames(GSE74639),
                                          rownames(GSE75367)))
length(intersect_genes)
saveRDS(intersect_genes,file="common_ctc_genes.rds")
write.csv(intersect_genes,file="common_ctc_genes.csv")
# Save only CTC Data
data <- cbind(  GSE51827[intersect_genes,],
                GSE55807[intersect_genes,],
                GSE60407[intersect_genes,],
                GSE67939[intersect_genes,1:15],
                GSE67980[intersect_genes,1:164],
                GSE74639[intersect_genes,],
                GSE75367[intersect_genes,])
new_col_names= c(
  c(paste("CTC_GSE51827", colnames(GSE51827) ,c(1:29),sep="_"))
  ,c(paste("CTC_GSE55807", colnames(GSE55807) ,c(1:6),sep="_"))
  ,c(paste("CTC_GSE60407", colnames(GSE60407) ,c(1:7),sep="_"))
  ,c(paste("CTC_GSE67939", colnames(GSE67939[,1:15]) ,c(1:15),sep="_"))
  ,c(paste("CTC_GSE67980", colnames(GSE67980[,1:164]) ,c(1:164),sep="_"))
  ,c(paste("CTC_GSE74639", colnames(GSE74639) ,c(1:16),sep="_"))
  ,c(paste("CTC_GSE75367", colnames(GSE75367) ,c(1:74),sep="_")))
colnames(data)<-new_col_names
dim(data)
#saveRDS(data,file='ctc__data.rds')


# Saving the Average Expression Per gene and Average Expression Per Sample
write.csv(colMeans(GSE51827[intersect_genes,]),file="GSE51827_sample_mean.csv")
write.csv(rowMeans(GSE51827[intersect_genes,]),file="GSE51827_gene_mean.csv")
write.csv(colMeans(GSE55807[intersect_genes,]),file="GSE55807_sample_mean.csv")
write.csv(rowMeans(GSE55807[intersect_genes,]),file="GSE55807_gene_mean.csv")
write.csv(colMeans(GSE60407[intersect_genes,]),file="GSE60407_sample_mean.csv")
write.csv(rowMeans(GSE60407[intersect_genes,]),file="GSE60407_gene_mean.csv")
write.csv(colMeans(GSE67939[intersect_genes,]),file="GSE67939_sample_mean.csv")
write.csv(rowMeans(GSE67939[intersect_genes,]),file="GSE67939_gene_mean.csv")
write.csv(colMeans(GSE67980[intersect_genes,]),file="GSE67980_sample_mean.csv")
write.csv(rowMeans(GSE67980[intersect_genes,]),file="GSE67980_gene_mean.csv")
write.csv(colMeans(GSE74639[intersect_genes,]),file="GSE74639_sample_mean.csv")
write.csv(rowMeans(GSE74639[intersect_genes,]),file="GSE74639_gene_mean.csv")
write.csv(colMeans(GSE75367[intersect_genes,]),file="GSE75367_sample_mean.csv")
write.csv(rowMeans(GSE75367[intersect_genes,]),file="GSE75367_gene_mean.csv")
write.csv(colMeans(pbmc_3k[intersect_genes,]),file="pbmc_3k_sample_mean.csv")
write.csv(rowMeans(pbmc_3k[intersect_genes,]),file="pbmc_3k_gene_mean.csv")
write.csv(colMeans(pbmc_6k[intersect_genes,]),file="pbmc_6k_sample_mean.csv")
write.csv(rowMeans(pbmc_6k[intersect_genes,]),file="pbmc_6k_gene_mean.csv")
write.csv(colMeans(GSE81861[intersect_genes,]),file="GSE81861_sample_mean.csv")
write.csv(rowMeans(GSE81861[intersect_genes,]),file="GSE81861_gene_mean.csv")
write.csv(Matrix::colMeans(EGAS00001002560[intersect_genes,]),file="EGAS00001002560_sample_mean.csv")
write.csv(Matrix::rowMeans(EGAS00001002560[intersect_genes,]),file="EGAS00001002560_gene_mean.csv")