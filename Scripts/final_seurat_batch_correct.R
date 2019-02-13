# ---
# Title: iCTC Project
# Description: Seurat Batch Processing
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
library(Seurat)

# Load the dataset
normalized_matrix_filter<-readRDS('../Data/orginal_normalized_ctc_blood_data.rds')
dim(normalized_matrix_filter)
# Seurat CCA Method Batch Effect
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
#CTC_GSE51827
GSE51827 <- CreateSeuratObject(raw.data = as.matrix(normalized_matrix_filter[,CTC_GSE51827_sample]))
GSE51827 <- NormalizeData(object=GSE51827)
GSE51827 <- FindVariableGenes(object=GSE51827, do.plot=FALSE)
GSE51827 <- ScaleData(object = GSE51827)
GSE51827@meta.data$tech <- "CTC_GSE51827"
#GSE51827 <- RenameCells(GSE51827, add.cell.id = "S1")
any(duplicated(x = colnames(x = GSE51827@data))) # check if there are duplicated cell names
dim(GSE51827@raw.data)
#GSE51827@var.genes <- genes

#CTC_GSE55807
GSE55807 <- CreateSeuratObject(raw.data = as.matrix(normalized_matrix_filter[,CTC_GSE55807_sample]))
GSE55807 <- NormalizeData(object=GSE55807)
GSE55807 <- FindVariableGenes(object=GSE55807, do.plot=FALSE)
GSE55807 <- ScaleData(object = GSE55807)
GSE55807@meta.data$tech <- "CTC_GSE55807"
#GSE55807 <- RenameCells(GSE55807, add.cell.id = "S2")
any(duplicated(x = colnames(x = GSE55807@data))) # check if there are duplicated cell names
dim(GSE55807@raw.data)
#GSE55807@var.genes <- genes

#CTC_GSE60407
GSE60407 <- CreateSeuratObject(raw.data = as.matrix(normalized_matrix_filter[,CTC_GSE60407_sample]))
GSE60407 <- NormalizeData(object=GSE60407)
GSE60407 <- FindVariableGenes(object=GSE60407, do.plot=FALSE)
GSE60407 <- ScaleData(object = GSE60407)
GSE60407@meta.data$tech <- "CTC_GSE60407"
#GSE60407 <- RenameCells(GSE60407, add.cell.id = "S3")
any(duplicated(x = colnames(x = GSE60407@data))) # check if there are duplicated cell names
dim(GSE60407@raw.data)
#GSE60407@var.genes <- genes

#CTC_GSE74639
GSE74639 <- CreateSeuratObject(raw.data = normalized_matrix_filter[,CTC_GSE74639_sample])
GSE74639 <- NormalizeData(object=GSE74639)
GSE74639 <- FindVariableGenes(object=GSE74639, do.plot=FALSE)
GSE74639 <- ScaleData(object = GSE74639)
GSE74639@meta.data$tech <- "CTC_GSE74639"
#GSE74639 <- RenameCells(GSE74639, add.cell.id = "S4")
any(duplicated(x = colnames(x = GSE74639@data))) # check if there are duplicated cell names
dim(GSE74639@raw.data)
#GSE74639@var.genes <- genes

#CTC_GSE67980
GSE67980 <- CreateSeuratObject(raw.data = normalized_matrix_filter[,CTC_GSE67980_sample])
GSE67980 <- NormalizeData(object=GSE67980)
GSE67980 <- FindVariableGenes(object=GSE67980, do.plot=FALSE)
GSE67980 <- ScaleData(object = GSE67980)
GSE67980@meta.data$tech <- "CTC_GSE67980"
#GSE67980 <- RenameCells(GSE67980, add.cell.id = "S5")
any(duplicated(x = colnames(x = GSE67980@data))) # check if there are duplicated cell names
dim(GSE67980@raw.data)
#GSE67980@var.genes <- genes

#CTC_GSE75367
GSE75367 <- CreateSeuratObject(raw.data = normalized_matrix_filter[,CTC_GSE75367_sample])
GSE75367 <- NormalizeData(object=GSE75367)
GSE75367 <- FindVariableGenes(object=GSE75367, do.plot=FALSE)
GSE75367 <- ScaleData(object = GSE75367)
GSE75367@meta.data$tech <- "CTC_GSE75367"
#GSE75367 <- RenameCells(GSE75367, add.cell.id = "S6")
any(duplicated(x = colnames(x = GSE75367@data))) # check if there are duplicated cell names
#GSE75367@var.genes <- genes
dim(GSE75367@raw.data)

#Blood_GSE81861
GSE81861 <- CreateSeuratObject(raw.data = normalized_matrix_filter[,Blood_GSE81861_sample])
GSE81861 <- NormalizeData(object=GSE81861)
#GSE81861 <- RenameCells(GSE81861, add.cell.id = "S7")
GSE81861 <- FindVariableGenes(object=GSE81861, do.plot=FALSE)
GSE81861 <- ScaleData(object = GSE81861)
GSE81861@meta.data$tech <- "Blood_GSE81861"
dim(GSE81861@raw.data)
any(duplicated(x = colnames(x = GSE81861@data))) # check if there are duplicated cell names
#GSE81861@var.genes <- genes

#Blood_Satijalab
Satijalab <- CreateSeuratObject(raw.data = normalized_matrix_filter[,Blood_Satijalab_sample])
Satijalab <- NormalizeData(object=Satijalab)
#Satijalab <- RenameCells(Satijalab, add.cell.id = "S8")
Satijalab <- FindVariableGenes(object=Satijalab, do.plot=FALSE)
Satijalab <- ScaleData(object = Satijalab)
Satijalab@meta.data$tech <- "Blood_Satijalab"
any(duplicated(x = colnames(x = Satijalab@data))) # check if there are duplicated cell names
#Satijalab@var.genes <- genes

#PBMC_5419
PBMC_5419 <- CreateSeuratObject(raw.data = normalized_matrix_filter[,Blood_PBMC_5419_sample])
PBMC_5419 <- NormalizeData(object=PBMC_5419)
#PBMC_5419 <- RenameCells(PBMC_5419, add.cell.id = "S9")
PBMC_5419 <- FindVariableGenes(object=PBMC_5419, do.plot=FALSE)
PBMC_5419 <- ScaleData(object = PBMC_5419)
PBMC_5419@meta.data$tech <- "Blood_PBMC_5419"
any(duplicated(x = colnames(x = PBMC_5419@data))) # check if there are duplicated cell names
#PBMC_5419@var.genes <- genes

#genomics_10x
genomics_10x <- CreateSeuratObject(raw.data = normalized_matrix_filter[,Blood_genomics_10x_mtx])
genomics_10x <- NormalizeData(object=genomics_10x)
#genomics_10x <- RenameCells(genomics_10x, add.cell.id = "S10")
genomics_10x <- FindVariableGenes(object=genomics_10x, do.plot=FALSE)
genomics_10x <- ScaleData(object = genomics_10x)
genomics_10x@meta.data$tech <- "Blood_genomics_10x"
any(duplicated(x = colnames(x = genomics_10x@data))) # check if there are duplicated cell names
#genomics_10x@var.genes <- genes

#Blood_GSE67980
Blood_GSE67980 <-CreateSeuratObject(raw.data = normalized_matrix_filter[,Blood_GSE67980_sample])
Blood_GSE67980 <- NormalizeData(object=Blood_GSE67980)
#Blood_GSE67980 <- RenameCells(Blood_GSE67980, add.cell.id = "S11")
Blood_GSE67980 <- FindVariableGenes(object=Blood_GSE67980, do.plot=FALSE)
Blood_GSE67980 <- ScaleData(object = Blood_GSE67980)
Blood_GSE67980@meta.data$tech <- "Blood_GSE67980"
any(duplicated(x = colnames(x = Blood_GSE67980@data))) # check if there are duplicated cell names
#Blood_GSE67980@var.genes <- genes

#CTC_GSE67939
GSE67939 <- CreateSeuratObject(raw.data = normalized_matrix_filter[,CTC_GSE67939_sample])
GSE67939 <- NormalizeData(object=GSE67939)
#GSE67939 <- RenameCells(GSE67939, add.cell.id = "S12")
GSE67939 <- FindVariableGenes(object=GSE67939, do.plot=FALSE)
GSE67939 <- ScaleData(object = GSE67939)
GSE67939@meta.data$tech <- "CTC_GSE67939"
any(duplicated(x = colnames(x = GSE67939@data))) # check if there are duplicated cell names


data1 <- list(GSE51827,
              GSE55807,
              GSE74639,
              GSE60407,
              GSE67939,
              GSE67980,
              GSE75367,
              GSE81861,
              Satijalab,
              PBMC_5419,
              genomics_10x,
              Blood_GSE67980)

genes.use <- c()
for (i in 1:length(data1)) {
  genes.use <- c(genes.use, head(rownames(data1[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(data1)) {
  genes.use <- genes.use[genes.use %in% rownames(data1[[i]]@scale.data)]
}
print(length(genes.use))
dim(normalized_matrix_filter)
# CCA
all_cca <- RunMultiCCA(object.list = data1,genes.use = genes.use, num.ccs = 2)
#Alignment
all_cca<- AlignSubspace(all_cca,reduction.type = "cca",grouping.var = "tech", dims.align = 1:2)

dim(all_cca@scale.data)