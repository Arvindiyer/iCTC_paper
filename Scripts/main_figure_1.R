# ---
# Title: iCTC Project
# Description: Main Figure-1
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/Data/')
#setwd('~/Data/')
set.seed(10)
library(readr)
library(ggpubr)
library(viridis)
library(superheat)
library(RColorBrewer)

normalized_data <- readRDS(file="ctc_blood_normalized_data.rds")

# Marker Genes
fibro_reduced_genes <- c("VIM","CD44","CD47","CD81","CKAP4" )
epi_reduced_genes <- c("KRT7","KRT8","KRT18","KRT19","CDH1","EPCAM","ERBB2","CAV2","DSP","TSPAN13","KLF6","AMACR","MUC1","CLDN4")
platelet_genes_reduced <- c("ITGB3","SELP","CD151","CD84","CD9","CD46","CD47","CD48","CD63")
t_b_cell_reduced <- c("CD27", "CD52","CMTM7","CD37","IL7R","CD27" )

# Subsmapling the blood data to sampe number of samples as CTC
dim(normalized_data)
length(grep("CTC_", colnames(normalized_data)))
length(grep("Blood_", colnames(normalized_data)))
temp<-311:13649
index<-sample(temp, 310)
ctc<-1:310
new_index<-c(ctc,index)
ctc_label<-rep('CTC',310)
blood_label<-rep('Blood',310)
conditon<-c(ctc_label,blood_label)
# Creating he data to plot 
new_mat_to_plot<-normalized_data[na.omit(match(epi_reduced_genes,rownames(normalized_data))),new_index]
new_mat_to_plot<-normalized_data[na.omit(match(fibro_reduced_genes,rownames(normalized_data))),new_index]
new_mat_to_plot<-normalized_data[na.omit(match(t_b_cell_reduced,rownames(normalized_data))),new_index]
new_mat_to_plot<-normalized_data[na.omit(match(platelet_genes_reduced,rownames(normalized_data))),new_index]

# Plot the heatmap Figure 1 A,B,C,D
superheat(log(new_mat_to_plot+1),
          left.label.col = "white",
          left.label.text.size = 3.5,
          bottom.label.col = "white",
          grid.hline = T,
          grid.vline = T,
          grid.hline.col = "black",
          grid.vline.col = "black",
          membership.cols = conditon,
          legend.width = 0.5,
          heat.lim = c(0,4),
          heat.pal = colorRampPalette(rev(brewer.pal(n = 5, name ="Spectral")))(100)
)

# Figure 1-G data extraction
# Load the ML used data
normalized_data <- readRDS('orginal_normalized_ctc_blood_data.rds')
ctc_markers<-c('KRT7','KRT8','KRT18','KRT19','CDH1','EPCAM','ERBB2','CDH11')
blood_data<-normalized_data[na.omit(match(ctc_markers,rownames(normalized_data))),311:13649]
dim(blood_data)
express=apply(blood_data,2,function(x) sum(x>0)>=1)
# Number of samples having atleast one ctc marker expressed.
sum(express)
percentage <- sum(express)*100/13339 
# Checkig the prediction for these samples
load('../../mnn_direct_predicted_84_new.Rdata')
rf_list <-seq(51, 84, by=7)
tblood_data<-list()
k=1
for (i in rf_list) {
  print(mnn_direct_predicted_84_new[[i]]$case)
  print(length(mnn_direct_predicted_84_new[[i]]$te_labels))
  tblood_data[[k]] <- mnn_direct_predicted_84_new[[i]]$te_labels
  k=k+1
}
t_blood<-unlist(tblood_data)
length(t_blood)
test<-t_blood[express]
# Prediction Summary
table(test)
#Figure 1-E
ml_plot_data <- read_csv("../../ml_ballon_plot_data.csv", 
                                col_types = cols(MNN = col_number(), 
                                                 Orginal = col_number(),
                                                 RCA = col_number(), 
                                                 Seurat = col_number()))

ml_plot_data<- as.data.frame(ml_plot_data)
rownames(ml_ballon_plot_data) <- ml_plot_data[, 1]
obj <- ml_plot_data[, -1]
pal <- viridisLite::viridis(6,option = "A",begin = .3,end = 1,direction = -1)
rownames(obj)<-c("CTC Breast #1 (Aceto N et al.)",
                 "CTC Breast #2 (Yu M et al.)",
                 "CTC Pancreas #3 (Ting DT et al.)",
                 "CTC Breast #4 (Sarioglu AF et al.)",
                 "CTC Prostrate #5 (Miyamoto DT et al.)",
                 "CTC Lung #6 (Zheng Y et al.)",
                 "CTC Breast #7 (Jordan NV et al.)",
                 "Blood #1 (Li H et al.)",
                 "Blood #2 (Grace XY Zheng et al.)",
                 "Blood #3 (Grace XY Zheng et al.)",
                 "Blood #4 (Miyamoto DT et al.)",
                 "Blood #5 (MGP van der Wijst, et al)")
superheat( obj,
           left.label.col = "white",
           left.label.text.col = "black",
           left.label.text.size = 3,
           bottom.label.text.size = 4,
           left.label.text.alignment = "left",
           bottom.label.col = "white",
           grid.hline = T,
           grid.vline = T,
           grid.hline.col = "white",
           grid.vline.col = "white",
           #legend.vspace = 0.1,
           left.label.size = 0.85,
           bottom.label.size = 0.1,
           heat.pal = pal,
           X.text = as.matrix(obj),
           X.text.size = 4,
           X.text.col = "black",
           legend = F)
#Figure 1-F
ml_plot_data <- read_csv("../../ml_box_plot_data.csv", 
                         col_types = cols(Accuracy = col_number(), 
                                          MCC = col_number(),
                                          Kappa = col_double(), 
                                          Method = col_character()))

ggboxplot(ml_plot_data,
          x = "Method",
          y = "Kappa",
          fill = "Method",
          bxp.errorbar = TRUE,
          notch = FALSE,
          palette = "Set2")+
  guides(size = FALSE)+
  theme(axis.text.y = element_text(family = "Helvitca",colour = "black",size = 12),
        axis.text.x = element_text(family = "Helvitca",colour = "black",size = 12,angle = 0, hjust = 0.5),
        legend.position = "None",
        legend.text = element_text(family = "Helvitca",color = "black",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"))+
  xlab('')+
  ylab('Cohen\'s Kappa Score')