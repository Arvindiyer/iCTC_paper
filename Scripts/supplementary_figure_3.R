# ---
# Title: iCTC Project
# Description: Supplementary Figure-3
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/Data/')
#setwd('~/Data/')
set.seed(10)

library(dplyr)
library(tidyr)
library(stringi)
library(stringr)
library(ggpubr)

#---
# Color Pallete  
# --
color_pal = c("#86BF4D", "#660000", "#B2DF8A", "#7570B3", "#A6761D", "#D95F02", "#E6AB02", "#E7298A", "#0020E9", "#770078", "#DB0000", "#800000")

#---
# Load Combined dataset  
# --
filter_normalized_data <- readRDS('../Data/orginal_normalized_ctc_blood_data.rds')

#---
# Scatter Plot : Orginal 
# --
tsne_data <- readRDS('../Data/final_orginal_tsne_data.rds')
plt_data = data.frame(  tsne1 = tsne_data$Y[,1],tsne2 = tsne_data$Y[,2], cluster=colnames(filter_normalized_data))
new_plt_data <- plt_data %>% separate(cluster, c("Type", "Study","Sample"), "_")
new_plt_data[new_plt_data=='genomics'] <- 'EGAS00001002560'
new_plt_data[new_plt_data=='Satijalab'] <- 'PBMC 3K'
new_plt_data[new_plt_data=='PBMC'] <- 'PBMC 6K'
ggscatter(x="tsne1",
          y="tsne2",
          data = new_plt_data,
          color = "Study",
          shape = "Type",
          palette = color_pal,
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

ctc_data <- new_plt_data[new_plt_data$Type =='CTC',]
ggscatter(x="tsne1",
          y="tsne2",
          data = ctc_data,
          color = "Study",
          shape = "Type",
          palette = "Set1",
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

blood_data <- new_plt_data[new_plt_data$Type =='Blood',]
ggscatter(x="tsne1",
          y="tsne2",
          data = blood_data,
          color = "Study",
          shape = "Type",
          palette = "Set1",
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

#---
# Scatter Plot : RCA 
# --
tsne_data <- readRDS('../Data/final_rca_tsne_data.rds')
plt_data = data.frame(  tsne1 = tsne_data$Y[,1],tsne2 = tsne_data$Y[,2], cluster=colnames(filter_normalized_data))
new_plt_data <- plt_data %>% separate(cluster, c("Type", "Study","Sample"), "_")
new_plt_data[new_plt_data=='genomics'] <- 'EGAS00001002560'
new_plt_data[new_plt_data=='Satijalab'] <- 'PBMC 3K'
new_plt_data[new_plt_data=='PBMC'] <- 'PBMC 6K'
ggscatter(x="tsne1",
          y="tsne2",
          data = new_plt_data,
          color = "Study",
          shape = "Type",
          palette = color_pal,
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

ctc_data <- new_plt_data[new_plt_data$Type =='CTC',]
ggscatter(x="tsne1",
          y="tsne2",
          data = ctc_data,
          color = "Study",
          shape = "Type",
          palette = "Set1",
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

blood_data <- new_plt_data[new_plt_data$Type =='Blood',]
ggscatter(x="tsne1",
          y="tsne2",
          data = blood_data,
          color = "Study",
          shape = "Type",
          palette = "Set1",
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))
#---
# Scatter Plot : MNN 
# --
tsne_data <- readRDS('../Data/final_mnn_tsne_data.rds')
plt_data = data.frame(  tsne1 = tsne_data$Y[,1],tsne2 = tsne_data$Y[,2], cluster=colnames(filter_normalized_data))
new_plt_data <- plt_data %>% separate(cluster, c("Type", "Study","Sample"), "_")
new_plt_data[new_plt_data=='genomics'] <- 'EGAS00001002560'
new_plt_data[new_plt_data=='Satijalab'] <- 'PBMC 3K'
new_plt_data[new_plt_data=='PBMC'] <- 'PBMC 6K'
ggscatter(x="tsne1",
          y="tsne2",
          data = new_plt_data,
          color = "Study",
          shape = "Type",
          palette = color_pal,
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))
ctc_data <- new_plt_data[new_plt_data$Type =='CTC',]
ggscatter(x="tsne1",
          y="tsne2",
          data = ctc_data,
          color = "Study",
          shape = "Type",
          palette = "Set1",
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

blood_data <- new_plt_data[new_plt_data$Type =='Blood',]
ggscatter(x="tsne1",
          y="tsne2",
          data = blood_data,
          color = "Study",
          shape = "Type",
          palette = "Set1",
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))


#---
# Scatter Plot : Seurat 
# --
tsne_data <- readRDS('../Data/final_seurat_tsne_data.rds')
plt_data = data.frame(  tsne1 = tsne_data$Y[,1],tsne2 = tsne_data$Y[,2], cluster=colnames(filter_normalized_data))
new_plt_data <- plt_data %>% separate(cluster, c("Type", "Study","Sample"), "_")
new_plt_data[new_plt_data=='genomics'] <- 'EGAS00001002560'
new_plt_data[new_plt_data=='Satijalab'] <- 'PBMC 3K'
new_plt_data[new_plt_data=='PBMC'] <- 'PBMC 6K'
ggscatter(x="tsne1",
          y="tsne2",
          data = new_plt_data,
          color = "Study",
          shape = "Type",
          palette = color_pal,
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

ctc_data <- new_plt_data[new_plt_data$Type =='CTC',]
ggscatter(x="tsne1",
          y="tsne2",
          data = ctc_data,
          color = "Study",
          shape = "Type",
          palette = "Set1",
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))

blood_data <- new_plt_data[new_plt_data$Type =='Blood',]
ggscatter(x="tsne1",
          y="tsne2",
          data = blood_data,
          color = "Study",
          shape = "Type",
          palette = "Set1",
          size = 1.2)+
  scale_shape_manual(values=c(20, 1))+
  scale_size_manual(values = c(1,5))+
  theme(axis.text = element_text(family = "Helvitca",size = 12),
        legend.position = "right",
        legend.text = element_text(family = "Helvitca",size = 12))+
  xlab('tSNE1')+
  ylab('tSNE2')+
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3)))
