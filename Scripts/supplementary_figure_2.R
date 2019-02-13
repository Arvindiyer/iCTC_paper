# ---
# Title: iCTC Project
# Description: Supplementary Figure-2
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/reboot/paper_work/Data/')
#setwd('~/Data/')
set.seed(10)

#---
# Average Gene Box Plot
# --

# Load the CTC Dataset
GSE51827<- read_csv(file = "GSE51827_gene_mean.csv")[,2]
GSE55807<- read_csv(file = "GSE55807_gene_mean.csv")[,2]
GSE60407<- read_csv(file = "GSE60407_gene_mean.csv")[,2]
GSE67939<- read_csv(file = "GSE67939_gene_mean.csv")[,2]
GSE67980<- read_csv(file = "GSE67980_gene_mean.csv")[,2]
GSE74639<- read_csv(file = "GSE74639_gene_mean.csv")[,2]
GSE75367<- read_csv(file = "GSE75367_gene_mean.csv")[,2]

study <- c(rep('GSE51827',15401),
           rep('GSE55807',15401),
           rep('GSE60407',15401),
           rep('GSE67939',15401),
           rep('GSE67980',15401),
           rep('GSE74639',15401),
           rep('GSE75367',15401))

average_read_per_gene <- c(GSE51827$x,
                           GSE55807$x,
                           GSE60407$x,
                           GSE67939$x,
                           GSE67980$x,
                           GSE74639$x,
                           GSE75367$x)

plot_data <- data.frame("average_read_per_gene" = log2(average_read_per_gene),
                           "Study" = study)

dim(plot_data)

ggboxplot(plot_data,
          x = "Study",
          y = "average_read_per_gene",
          fill = "Study",
          bxp.errorbar = TRUE,
          notch = FALSE,
          palette = "Set2")+
  guides(size = FALSE)+
  theme(plot.title = element_text(color="black", size=14, hjust = 0.5),
        axis.text.y = element_text(family = "Helvitca",colour = "black",size = 10),
        axis.text.x = element_text(family = "Helvitca",colour = "black",size = 10,angle = 0, hjust = 0.5),
        legend.position = "None",
        legend.text = element_text(family = "Helvitca",color = "black",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"))+
  xlab('')+
  ylab('log2 average reads per gene')+
  ggtitle("CTC studies")

# Load the Blood Dataset
GSE67980<- read_csv(file = "GSE67980_gene_mean.csv")[,2]
pbmc_3k<- read_csv(file="pbmc_3k_gene_mean.csv")[,2]
GSE81861<- read_csv(file="GSE81861_gene_mean.csv")[,2]
pbmc_6k<- read_csv(file="pbmc_6k_gene_mean.csv")[,2]
EGAS00001002560<- read_csv(file="EGAS00001002560_gene_mean.csv")[,2]

study <- c(rep('GSE67980',15401),
           rep('PBMC 3K',15401),
           rep('GSE81861',15401),
           rep('PBMC 6K',15401),
           rep('EGAS00001002560',15401))

average_read_per_gene <- c(GSE67980$x,
                           pbmc_3k$x,
                           GSE81861$x,
                           pbmc_6k$x,
                           EGAS00001002560$x)

plot_data <- data.frame("average_read_per_gene" = log2(average_read_per_gene),
                        "Study" = study)

dim(plot_data)

ggboxplot(plot_data,
          x = "Study",
          y = "average_read_per_gene",
          fill = "Study",
          bxp.errorbar = TRUE,
          notch = FALSE,
          palette = "Set2")+
  guides(size = FALSE)+
  theme(plot.title = element_text(color="black", size=14, hjust = 0.5),
        axis.text.y = element_text(family = "Helvitca",colour = "black",size = 10),
        axis.text.x = element_text(family = "Helvitca",colour = "black",size = 10,angle = 0, hjust = 0.5),
        legend.position = "None",
        legend.text = element_text(family = "Helvitca",color = "black",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"))+
  xlab('')+
  ylab('log2 average reads per gene')+
  ggtitle("Blood studies")


#---
# Average Sample Box Plot
# --

# Load the CTC Dataset
GSE51827<- read_csv(file = "GSE51827_sample_mean.csv")[,2]
GSE55807<- read_csv(file = "GSE55807_sample_mean.csv")[,2]
GSE60407<- read_csv(file = "GSE60407_sample_mean.csv")[,2]
GSE67939<- read_csv(file = "GSE67939_sample_mean.csv")[,2]
GSE67980<- read_csv(file = "GSE67980_sample_mean.csv")[,2]
GSE74639<- read_csv(file = "GSE74639_sample_mean.csv")[,2]
GSE75367<- read_csv(file = "GSE75367_sample_mean.csv")[,2]

study <- c(rep('GSE51827',29),
           rep('GSE55807',6),
           rep('GSE60407',7),
           rep('GSE67939',17),
           rep('GSE67980',169),
           rep('GSE74639',16),
           rep('GSE75367',74))

average_read_per_sample <- c(GSE51827$x,
                           GSE55807$x,
                           GSE60407$x,
                           GSE67939$x,
                           GSE67980$x,
                           GSE74639$x,
                           GSE75367$x)

plot_data <- data.frame("average_read_per_sample" = log2(average_read_per_sample+1),
                        "Study" = study)

dim(plot_data)

ggboxplot(plot_data,
          x = "Study",
          y = "average_read_per_sample",
          fill = "Study",
          bxp.errorbar = TRUE,
          notch = FALSE,
          palette = "Set2")+
  guides(size = FALSE)+
  theme(plot.title = element_text(color="black", size=14, hjust = 0.5),
        axis.text.y = element_text(family = "Helvitca",colour = "black",size = 10),
        axis.text.x = element_text(family = "Helvitca",colour = "black",size = 10,angle = 0, hjust = 0.5),
        legend.position = "None",
        legend.text = element_text(family = "Helvitca",color = "black",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"))+
  xlab('')+
  ylab('log2 average reads per sample')+
  ggtitle("CTC studies")

# Load the Blood Dataset
GSE67980<- read_csv(file = "GSE67980_sample_mean.csv")[,2]
pbmc_3k<- read_csv(file="pbmc_3k_sample_mean.csv")[,2]
GSE81861<- read_csv(file="GSE81861_sample_mean.csv")[,2]
pbmc_6k<- read_csv(file="pbmc_6k_sample_mean.csv")[,2]
EGAS00001002560<- read_csv(file="EGAS00001002560_sample_mean.csv")[,2]

study <- c(rep('GSE67980',169),
           rep('PBMC 3K',142),
           rep('GSE81861',2700),
           rep('PBMC 6K',5419),
           rep('EGAS00001002560',28855))

average_read_per_gene <- c(GSE67980$x,
                           pbmc_3k$x,
                           GSE81861$x,
                           pbmc_6k$x,
                           EGAS00001002560$x)

plot_data <- data.frame("average_read_per_gene" = log2(average_read_per_gene),
                        "Study" = study)

dim(plot_data)

ggboxplot(plot_data,
          x = "Study",
          y = "average_read_per_gene",
          fill = "Study",
          bxp.errorbar = TRUE,
          notch = FALSE,
          palette = "Set2")+
  guides(size = FALSE)+
  theme(plot.title = element_text(color="black", size=14, hjust = 0.5),
        axis.text.y = element_text(family = "Helvitca",colour = "black",size = 10),
        axis.text.x = element_text(family = "Helvitca",colour = "black",size = 10,angle = 0, hjust = 0.5),
        legend.position = "None",
        legend.text = element_text(family = "Helvitca",color = "black",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"))+
  xlab('')+
  ylab('log2 average reads per sample')+
  ggtitle("Blood studies")