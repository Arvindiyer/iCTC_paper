# ---
# Title: iCTC Project
# Description: E:M analysis script
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# version: 0.2
# Updates:
# 1] Removed overlapping study data
# 2] Chnaged the cutoff and no-need to resample the blood data
# 3] Solved few bugs and refactor the code
# ---

# Set working Directory and a seed
setwd('~/final_github_iCTC/')
#setwd('~/Data/')
set.seed(10)
# load utility function
source('scripts/functions.R')

library(readr)
library(ggpubr)
library(viridis)
library(superheat)
library(RColorBrewer)
library(pheatmap)

# Load ctc data only
data <- readRDS('dataset/ctc_raw_count_data.rds')
dim(data)
raw.counts<-apply(data,2,function(x) {storage.mode(x) <- 'integer'; x})
dim(raw.counts)
# Grep the total CTC and removed overlapping study data
length(grep("CTC_", colnames(raw.counts)))
CTC_GSE51827_sample=grep("CTC_GSE51827", colnames(data))
CTC_GSE55807_sample=grep("CTC_GSE55807", colnames(data))
CTC_GSE60407_sample=grep("CTC_GSE60407", colnames(data))
CTC_GSE67939_sample=grep("CTC_GSE67939", colnames(data))
CTC_GSE67980_sample=grep("CTC_GSE67980", colnames(data))
CTC_GSE74639_sample=grep("CTC_GSE74639", colnames(data))
CTC_GSE75367_sample=grep("CTC_GSE75367", colnames(data))
CTC_GSE109761_sample = grep("CTC_GSE109761",colnames(data))
CTC_GSE86978_sample = grep("CTC_GSE86978",colnames(data))
CTC_GSE38495_sample = grep("CTC_GSE38495",colnames(data))
index <- index <- c(CTC_GSE51827_sample,
                    CTC_GSE55807_sample,
                    CTC_GSE60407_sample,
                    CTC_GSE67939_sample,
                    CTC_GSE67980_sample,
                    CTC_GSE74639_sample,
                    CTC_GSE75367_sample,
                    CTC_GSE109761_sample,
                    CTC_GSE86978_sample,
                    CTC_GSE38495_sample)

reduced_data <- raw.counts[,index]
dim(reduced_data)
rm(data)
# Apply cell filtering to the dataset
cell_filter_data <- cell_filter(reduced_data,5)
dim(cell_filter_data)
# Apply gene filtering to the dataset
gene_filter_data <- gene_filter(cell_filter_data,5,10)
dim(gene_filter_data)
# Load Epithelial Genes
epi_genes <- readRDS('dataset/final_epitheilial_genes.rds')
index<- na.omit(match(toupper(epi_genes),rownames(gene_filter_data)))
new_genes<-na.omit(rownames(gene_filter_data)[index])
epi_genes <- new_genes
mesn_genes <- readRDS('dataset/final_mesenchymal_genes.rds')
index<- na.omit(match(toupper(mesn_genes),rownames(gene_filter_data)))
new_genes<-na.omit(rownames(gene_filter_data)[index])
mesn_genes <- new_genes
genes<-read.csv('dataset/Cancer_stem_marker_genes.csv')
index<- na.omit(match(toupper(genes$genes),rownames(gene_filter_data)))
new_genes<-na.omit(rownames(gene_filter_data)[index])
csc_genes <- new_genes
length(unique(epi_genes))
length(unique(mesn_genes))
length(unique(csc_genes))
common1= intersect(unique(epi_genes),unique(csc_genes))
common2= intersect(unique(mesn_genes),unique(csc_genes))
common3= intersect(unique(mesn_genes),unique(epi_genes))
match1=match(common1,unique(csc_genes))
match2=match(common2,unique(csc_genes))
common1= intersect(unique(epi_genes),unique(csc_genes))
common2= intersect(unique(mesn_genes),unique(csc_genes))
common3= intersect(unique(mesn_genes),unique(epi_genes))
match1=match(common1,unique(csc_genes))
match2=match(common2,unique(csc_genes))
csc_genes_reduced <- unique(csc_genes[-match1])
csc_genes_reduced <- unique(csc_genes_reduced[-match2])
length(unique(epi_genes))
length(unique(mesn_genes))
length(unique(csc_genes_reduced))
total<-c(epi_genes,mesn_genes,csc_genes_reduced)
length(total)
length(unique(total))
gene_filter_data_reduced <- gene_filter_data[unique(total),]
dim(gene_filter_data_reduced)
cell_filter_data <- cell_filter(gene_filter_data_reduced,percentage = 10)
dim(cell_filter_data)
# Apply gene filtering to the dataset (30% of 550)
gene_filter_data <- gene_filter(cell_filter_data,5,165)
dim(gene_filter_data)
normalized_data <- normalization(gene_filter_data,method = "median")
new_mat_to_plot_gene<-unique(normalized_data)
saveRDS(normalized_data,file = "dataset/final_marker_normalized_matrix_v2.rds")

length(na.omit(match(unique(epi_genes),rownames(new_mat_to_plot_gene))))
length(na.omit(match(unique(mesn_genes),rownames(new_mat_to_plot_gene))))
length(na.omit(match(unique(csc_genes_reduced),rownames(new_mat_to_plot_gene))))
normalized_matrix_filter<-normalized_data
#normalized_matrix_filter <- readRDS("dataset/final_marker_normalized_matrix.rds")
CTC_GSE51827_sample= length(grep("CTC_GSE51827", colnames(normalized_matrix_filter)))
CTC_GSE55807_sample= length(grep("CTC_GSE55807", colnames(normalized_matrix_filter)))
CTC_GSE60407_sample= length(grep("CTC_GSE60407", colnames(normalized_matrix_filter)))
CTC_GSE67939_sample= length(grep("CTC_GSE67939", colnames(normalized_matrix_filter)))
CTC_GSE67980_sample= length(grep("CTC_GSE67980", colnames(normalized_matrix_filter)))
CTC_GSE74639_sample= length(grep("CTC_GSE74639", colnames(normalized_matrix_filter)))
CTC_GSE75367_sample= length(grep("CTC_GSE75367", colnames(normalized_matrix_filter)))
CTC_GSE109761_sample= length(grep("CTC_GSE109761",colnames(normalized_matrix_filter)))
CTC_GSE86978_sample = length(grep("CTC_GSE86978",colnames(normalized_matrix_filter)))
CTC_GSE38495_sample = length(grep("CTC_GSE38495",colnames(normalized_matrix_filter)))

col_label<-c(rep("CTC Breast #1 (Aceto N et al.)",CTC_GSE51827_sample),
             rep("CTC Breast #2 (Yu M et al.)",CTC_GSE55807_sample),
             rep("CTC Pancreas #3 (Ting DT et al.)",CTC_GSE60407_sample),
             rep("CTC Breast #4 (Sarioglu AF et al.)",CTC_GSE67939_sample),
             rep("CTC Prostrate #5 (Miyamoto DT et al.)",CTC_GSE67980_sample),
             rep("CTC Lung #6 (Zheng Y et al.)",CTC_GSE74639_sample),
             rep("CTC Breast #7 (Jordan NV et al.)",CTC_GSE75367_sample),
             rep("CTC Breast #8 (Barbara MS et al.)",CTC_GSE109761_sample),
             rep("CTC Breast #9 (Aceto N et al)",CTC_GSE86978_sample),
             rep("CTC Melanoma #10 (Daniel R et al)",CTC_GSE38495_sample))

dim(normalized_matrix_filter)
row_annotation = c(rep('Epithelial',16),rep('Mesenchymal',39),rep('Cancer Stem Cell',26))
length(row_annotation)
length(col_label)
filter_epi_genes <- rownames(normalized_data)[1:16]
#save(filter_epi_genes,file = "dataset/new_final_filter_epitheilal_genes.Rdata")
saveRDS(filter_epi_genes,file = "dataset/new_final_filter_epitheilal_genes.rds")
#write.csv(filter_epi_genes,file="dataset/new_final_filter_epitheilal_genes.csv")
filter_mesn_genes <- rownames(normalized_data)[17:55]
#save(filter_mesn_genes,file = "dataset/new_final_filter_mesenchymal_genes.Rdata")
saveRDS(filter_mesn_genes,file = "dataset/new_final_filter_mesenchymal_genes.rds")
#write.csv(filter_mesn_genes,file="dataset/new_final_filter_mesenchymal_genes.csv")
filter_csc_genes <- rownames(normalized_data)[56:81]
#save(filter_csc_genes,file = "dataset/new_final_filter_cancer_stem_genes.Rdata")
saveRDS(filter_csc_genes,file = "dataset/new_final_filter_cancer_stem_genes.rds")
#write.csv(filter_csc_genes,file="dataset/new_final_filter_cancer_stem_genes.csv")
# Calculating the z_score
z_score <-(get_Zscore((log(normalized_data+1))))
dim(na.omit(z_score))
annotation_col <- data.frame(Study = factor(col_label))
rownames(annotation_col) <- colnames(z_score)
dim(annotation_col)
annotation_rows <- data.frame(Marker = factor(row_annotation))
rownames(annotation_rows)<-rownames(normalized_matrix_filter)
dim(annotation_rows)
test1<-get_stouffer(z_score[1:16,],2)
test2<-get_stouffer(z_score[17:55,],2)
test3<-get_stouffer(z_score[56:81,],2)
new_data <- data.frame('Cancer'=test3,'Epitheial'=test1,'Mesenchymal'=test2)
dim(new_data)
dim(z_score)
cor(new_data$Epitheial,new_data$Mesenchymal)
cor(new_data$Mesenchymal,new_data$Cancer)
cor(new_data$Epitheial,new_data$Cancer)
tol11qualitative=c("#332288",
                   "#6699CC",
                   "#88CCEE",
                   "#44AA99",
                   "#117733",
                   "#999933",
                   "#DDCC77",
                   "#661100",
                   "#CC6677",
                   "#882255",
                   "#AA4499")
ann_colors = list(Study = c("CTC Breast #1 (Aceto N et al.)"="#332288",
                            "CTC Breast #2 (Yu M et al.)"="#6699CC",
                            "CTC Pancreas #3 (Ting DT et al.)"="#88CCEE",
                            "CTC Breast #4 (Sarioglu AF et al.)"="#44AA99",
                            "CTC Prostrate #5 (Miyamoto DT et al.)"="#117733",
                            "CTC Lung #6 (Zheng Y et al.)"="#999933",
                            "CTC Breast #7 (Jordan NV et al.)"="#DDCC77",
                            "CTC Breast #8 (Barbara MS et al.)"="#661100",
                            "CTC Breast #9 (Aceto N et al)"="#882255",
                            "CTC Melanoma #10 (Daniel R et al)"="#AA4499"))
matrix_data<-new_data
dim(matrix_data)
head(matrix_data)
matrix_data$Study <- col_label
ggscatter(round_df(matrix_data,2),
          "Epitheial",
          "Mesenchymal",
          color = "Study",
          shape = "Study",
          #          palette = tol11qualitative,
          size = 1.5)+
  #  scale_color_jcolors(palette = "pal8")+
  scale_shape_manual(values=c(11:22))+
  scale_color_manual(values = ann_colors$Study)+
  theme(axis.text = element_text(family = "Helvitca",size = 10),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvitca",size = 8))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  annotate("text", x = 3, y = 1.5, label = "bold(italic(r)) == -0.937",parse = TRUE)+
  xlab('Epithelial Singature')+
  ylab('Mesenchymal Singature')
ggplot(round_df(matrix_data,2), aes(x=Epitheial, y=Mesenchymal,color=Study,shape=Study)) +
  facet_wrap(. ~ Study) +
  scale_color_manual(values = ann_colors$Study)+
  #  scale_color_manual(values =tol11qualitative)+
  scale_shape_manual(values=c(11:22))+
  geom_point() +
  theme_bw()+
  theme(axis.text = element_text(family = "Helvitca",size = 10),
        legend.position = "right",
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        #        strip.text = element_text(family = "Helvitca",size = 5),
        legend.text = element_text(family = "Helvitca",size = 8))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  xlab('Epithelial Singature')+
  ylab('Mesenchymal Singature')

# Ordering based on E:M
E<- matrix_data$Epitheial+abs(min(matrix_data$Epitheial))+1
M<- matrix_data$Mesenchymal+abs(min(matrix_data$Mesenchymal))+1
E_M <- E/M
matrix_data$E_M <-E_M
dim(matrix_data)
order_matrix_data <- matrix_data[order(-matrix_data$E_M),]
dim(order_matrix_data)
head(order_matrix_data)
save(order_matrix_data,file = "dataset/new_final_ctc_signature_E_M_ordering.Rdata")
saveRDS(order_matrix_data,file = "dataset/new_final_ctc_signature_E_M_ordering.rds")
dim(new_mat_to_plot_gene)
EM_order_data <- new_mat_to_plot_gene[1:59,order(-matrix_data$E_M)]
dim(EM_order_data)
save(EM_order_data,file = "dataset/new_final_singature_epi_mesn_expression.Rdata")
saveRDS(EM_order_data,file = "dataset/new_final_singature_epi_mesn_expression.rds")
dim(EM_order_data)
data <- EM_order_data
for (i in c(1:5))
{
  data <- movingAverageByCol(t(EM_order_data),10)
}
dim(data)
annotation_col <- data.frame(Study = factor(order_matrix_data$Study))
rownames(annotation_col) <- colnames(EM_order_data)
pheatmap(log(t(data)+1),
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_cols = F,
         fontsize_row = 6,
         fontsize = 8,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_row = annotation_rows,
         cluster_rows = F,
         #         cutree_rows = F,
         scale = "none",
         #         treeheight_row = 0,
         color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100))
