# ---
# Title: iCTC Project
# Description: Survival analysis of Lung Data
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/Data/')
#setwd('~/Data/')
set.seed(10)

# Library
library(survival)
library(survminer)
library(survcomp)
library(randomForestSRC)

#---
# Utility Functions
#---

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
# Correlation function which return the correlation matrix #
get_corr <- function(data1,data2){
  data1[data1 == 0]<-1
  data2[data2 == 0]<-1
  correlation = cor(data1,data2,method = "pearson")
  return (correlation)
}

# Z-score compuation which return the z-score matrix #  
get_Zscore <- function(data){
  mean=colMeans(data)
  std=apply(data,2,sd)
  #return (t((t(data1)-mean)/std))
  return(scale(data,scale=std,center=mean))
}
# Stouffer score compuation which return the score matrix #    
get_stouffer <- function(data1,axis){
  return(apply(data1,axis,function(x){sum(x)/sqrt(length(x))}))
}

#---
# Utility Functions Ends
#---

# # Load the Data
# load('updated_mnn_only_ctc_expression_data.Rdata')
# filter_data<- final_run_mnn
# dim(filter_data)
# # 
# # # Load lung Data
# load('lung_final_mRNA_raw_count.data')
# dim(lung_mRNA)
# # 
# # # Making the CTC same that of lung
# ctc_mRNA <- filter_data[match(rownames(lung_mRNA),rownames(filter_data)),]
# dim(ctc_mRNA)
# save(ctc_mRNA,file = "ctc_lung_mRNA.Rdata")
# # 
# # # Calculating the stouffer score for lung samples
# # # Normalizing the lung Data
# lung_mRNA_normalized <- normalization(lung_mRNA,method = "median")
# dim(lung_mRNA_normalized)
# dim(ctc_mRNA)
# tissuse_corr = get_corr(lung_mRNA_normalized,ctc_mRNA)
# dim(tissuse_corr)
# tissuse_zscore = get_Zscore(tissuse_corr)
# dim(tissuse_zscore)
# tissuse_stouffer=get_stouffer(tissuse_zscore,1)
# length(tissuse_stouffer)
# # #Store the stouffer score
# lung_stouffer <- data.frame('feature'=colnames(lung_mRNA_normalized),'stouffer'=tissuse_stouffer)
# save(lung_stouffer,file = "lung_stouffer_mnn.Rdata")
# 
# #Creating Bootsrap 100 interation samples
# bootstrap_train<-list()
# bootstrap_test<-list()
# for (i in c(1:100)) {
#  features<-unique(colnames(lung_mRNA_normalized))
#  test<-sample(colnames(lung_mRNA_normalized),97)
#  train<-setdiff(features,test)
#  bootstrap_train[i]=list(train)
#  bootstrap_test[i]=list(test)
# }
# save(bootstrap_train,file='lung_train_samples.RData')
# save(bootstrap_test,file='lung_test_samples.RData')



# Actual Survival Analysis
#Clear workspace
#rm(list=ls(all=TRUE))
# Load the lung mRNA Data
load('lung_final_mRNA_raw_count.data')
duplicated.columns <- duplicated(t(lung_mRNA))
lung_mRNA <- lung_mRNA[, !duplicated.columns]
dim(lung_mRNA)
mRNA <- normalization(lung_mRNA,method = "median")
mRNA<-t(mRNA)
dim(mRNA)
# Load Clincial Data
load('lung_final_clincial.data')
#dim(unique(lung_clincial))
lung_clincial<-as.data.frame(unique(lung_clincial))
#lung_clincial$gender <- factor(lung_clincial$gender)
lung_clincial[lung_clincial$gender=='male',]$gender<-0
lung_clincial[lung_clincial$gender=='female',]$gender<-1
lung_clincial$gender <- as.numeric(as.character(lung_clincial$gender))
lung_clincial$stage<-as.factor(lung_clincial$stage)
str(lung_clincial)

# Load Survival
load('lung_final_survival.data')
#dim(unique(lung_clincial))
lung_survival<-unique(lung_survival)
str(lung_survival)
#lung_survival$status<-factor(lung_survival$status)
lung_survival[lung_survival$status=='alive',]$status<-0
lung_survival[lung_survival$status=='dead',]$status<-1
lung_survival$time<-lung_survival$time
#lung_survival[lung_survival$time<0,]$time=0
lung_survival$status <- as.numeric(as.character(lung_survival$status))
str(lung_survival)

# Load Stouffer data
load('lung_stouffer_mnn.Rdata')
lung_stouffer<-unique(lung_stouffer)
dim(lung_stouffer)

# Clincal+Stouffer
lung_ctc_clincal<-data.frame(feature=lung_clincial$feature,stouffer=lung_stouffer$stouffer,age=lung_clincial$age,gender=lung_clincial$gender,stage=lung_clincial$stage)
#lung_ctc_clincal$stage<-as.factor(lung_ctc_clincal$stage)
str(lung_ctc_clincal)


#mRNA+Stouffer
dim(mRNA)
dim(lung_stouffer)
mRNA_stouffer <-as.data.frame(mRNA)
mRNA_stouffer$stouffer <- lung_stouffer$stouffer
dim(mRNA_stouffer)

#mRNA+Stouffer+clinical
mRNA_clincial_stouffer<- mRNA_stouffer
mRNA_clincial_stouffer$age <- lung_ctc_clincal$age
mRNA_clincial_stouffer$gender <- lung_ctc_clincal$gender
#mRNA_clincial_stouffer$score <- lung_ctc_clincal$score
mRNA_clincial_stouffer$stage <- lung_ctc_clincal$stage
dim(mRNA_clincial_stouffer)
# Load Bootstrap test and train
load('lung_train_samples.RData')
load('lung_test_samples.RData')

c_index = list()
for (i in c(1:length(bootstrap_train))) {
  train_index = match(bootstrap_train[i][[1]],lung_clincial$feature)
  test_index = match(bootstrap_test[i][[1]],lung_clincial$feature)
  #print(test_index)  
  # survival time and status data
  #trainData = as.data.frame(unclass(lung_ctc_clincal[train_index,2:5]))
  #trainData = as.data.frame(mRNA[train_index,1:ncol(mRNA)])
  #trainData = as.data.frame(data.frame(mRNA_stouffer[train_index,1:ncol(mRNA_stouffer)]))
  trainData = as.data.frame(mRNA_clincial_stouffer[train_index,1:ncol(mRNA_clincial_stouffer)])
  #trainData = as.data.frame(unclass(mRNA_clincial_stouffer[train_index,1:ncol(mRNA_clincial_stouffer)]))
  #trainData = lung_stouffer[train_index,]$stouffer
  trainSurvStatus = lung_survival$status[train_index]
  trainSurvTime = lung_survival$time[train_index]
  
  # Creating survival train data frame
  #surv_data = as.data.frame(unclass(data.frame(time=trainSurvTime,status=trainSurvStatus, stouffer=trainData)))
  surv_data = as.data.frame(data.frame(time=trainSurvTime,status=trainSurvStatus, trainData))
  
  # # Creating survival test data frame
  #testData = as.data.frame(mRNA[test_index,1:ncol(mRNA)])
  #testData = as.data.frame(data.frame(mRNA_stouffer[test_index,1:ncol(mRNA_stouffer)]))
  testData = as.data.frame(data.frame(mRNA_clincial_stouffer[test_index,1:ncol(mRNA_clincial_stouffer)]))
  #testData = as.data.frame(unclass(data.frame(mRNA_clincial_stouffer[test_index,1:ncol(mRNA_clincial_stouffer)])))
  #testData = as.data.frame(unclass(lung_ctc_clincal[test_index,2:5]))
  #testData = lung_stouffer[test_index,]$stouffer
  #testData = data.frame(stouffer=testData)
  #testData = as.data.frame(lung_ctc_clincal[test_index,2:5])
  testSurvStatus = lung_survival$status[test_index]
  testSurvTime = lung_survival$time[test_index]
  #levels(testData$gender)=levels(surv_data$gender)
  rfsrc.model.fit <- rfsrc(Surv(time,status) ~ ., data = surv_data, ntree=1000, na.action="na.impute", splitrule="logrank", nsplit=1, importance="random", seed=-1)
  predictedResponse <- predict.rfsrc(rfsrc.model.fit, newdata = testData, na.action="na.impute")$predicted
  concordance <- concordance.index(predictedResponse, testSurvTime, testSurvStatus)$c.index
  # cox.model.fit<-coxph(Surv(time,status) ~ ., data = surv_data)
  # concordance <- unlist(survConcordance(Surv(testSurvTime, testSurvStatus) ~ predict(cox.model.fit, testData), testData)$concordance)
  print(i)
  #Appending c_index value
  c_index<-c(c_index,concordance)
}
#save(c_index,file = "lung_stouffer_c_index_updated.Rdata")
#save(c_index,file = "lung_stouffer_c_index_updated_mnn.Rdata")
#save(c_index,file = "lung_clinical_c_index_updated.Rdata")
#save(c_index,file = "lung_clinical_stouffer_c_index_updated.Rdata")
#save(c_index,file = "lung_clinical_stouffer_c_index_updated_mnn.Rdata")
#save(c_index,file = "lung_mRNA_index_updated.Rdata")
#save(c_index,file = "lung_mRNA_stouffer_index_updated.Rdata")
#save(c_index,file = "lung_mRNA_stouffer_index_updated_mnn.Rdata")
#save(c_index,file = "lung_mRNA_clinical_stouffer_index_updated.Rdata")
save(c_index,file = "lung_mRNA_clinical_stouffer_index_updated_mnn.Rdata")

#Creating a combined lung cindex object
load("lung_clinical_c_index_updated.Rdata")
clinical<-c_index
load("lung_stouffer_c_index_updated_mnn.Rdata")
ctc_stouffer <- c_index
load("lung_mRNA_stouffer_index_updated_mnn.Rdata")
ctc_mRNA <- c_index
load("lung_clinical_stouffer_c_index_updated_mnn.Rdata")
ctc_clinical <- c_index
load("lung_mRNA_stouffer_index.Rdata")
mRNA <- c_index
load("lung_mRNA_clinical_stouffer_index_updated_mnn.Rdata")
ctc_mRNA_clinical <- c_index
# # Combine
c_index_data<-c(unlist(ctc_stouffer),unlist(clinical),unlist(mRNA),unlist(ctc_clinical),unlist(ctc_mRNA),unlist(ctc_mRNA_clinical))
label<-c(rep('Stouffer',100),rep('Clinical',100),rep('mRNA',100),rep('Clinical+',100),rep('mRNA+',100),rep('All',100))
length(c_index_data)
length(label)
lung_c_index_data<-data.frame('c-index values'=c_index_data,'methods'=label)
save(lung_c_index_data,file = "updated_lung_c_index_data_mnn.Rdata")
lung_Cummilative <- data.frame('CTC'=unlist(ctc_stouffer),'Clinical'=unlist(clinical),'mRNA'=unlist(mRNA),'CTC_Clinical'=unlist(ctc_clinical),'CTC_mRNA'=unlist(ctc_mRNA),'CTC_mRNA_Clinical'=unlist(ctc_mRNA_clinical))
save(lung_Cummilative,file = "updated_lung_box_plot_mnn.Rdata")
write.csv(lung_Cummilative,file = "updated_lung_box_plot_mnn.csv")
summary(lung_Cummilative)