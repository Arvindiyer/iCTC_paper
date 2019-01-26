# ---
# Title: iCTC Project
# Description: Raw Data and Preprocessing
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
#setwd('~/iCTC/reboot/paper_work/Data/')
setwd('~/Data/')
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
saveRDS(intersect_genes,file="common_ctc_blood_genes.rds")

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