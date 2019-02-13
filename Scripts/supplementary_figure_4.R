# ---
# Title: iCTC Project
# Description: Supplementary Figure-4
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/reboot/paper_work/Data/')
#setwd('~/Data/')
set.seed(10)

library(pheatmap)
library(readr)
library(ggpubr)
library(viridis)

# Load the ML prediction heat map data
load('cor_prediction_results.Rdata')
pheatmap(cor_prediction_results,
         cluster_row = T,
         cluster_cols = T,
         treeheight_row = 0,
         treeheight_col = 0)
# Load the Kappa Score Data 
ml_plot_data <- read_csv("ml_supplementary_data.csv")
head(ml_plot_data)
ggboxplot(ml_plot_data,
          x = "Method",
          y = "Kappa_Test",
          facet.by = "`Data Type`",
          fill = "Method",
          bxp.errorbar = TRUE,
          notch = FALSE,
          palette = "Set1")+
  guides(size = FALSE)+
  theme(axis.text.y = element_text(family = "Helvitca",colour = "black",size = 12),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(family = "Helvitca",color = "black",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"))+
  xlab('')+
  ylab('Cohen\'s Kappa Score')