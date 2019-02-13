# ---
# Title: iCTC Project
# Description: TCGA Data Processing
# Authors: Arvind Iyer <arvind16122@iiitd.ac.in>, Krishan Gupta <krishang@iitd.ac.in>, Shreya Sharma <shreya15096@iiitd.ac.in> 
# Corresponding Author: <Debarka Sengputa<debarka@iiitd.ac.in>
# Feel free to get in touch with us as we would love to talk and discuss science :):)
# ---

# Set working Directory and a seed
setwd('~/iCTC/Data/')
#setwd('~/Data/')
set.seed(10)
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
#---
# Utility Function End
#---

library(recount)

# Read only CTC Data genes
load('../../only_ctc_updated_filter.Rdata')
dim(filter_data)
genes<-rownames(filter_data)
#---
# GBM Data preprocessing
#---

#load GBM dataset 
load('../../rse_gene_brain.Rdata')
dim(rse_gene)
rse<-scale_counts(rse_gene)
rm(rse_gene)
dim(rse)
table(rse@colData@listData$gdc_cases.project.project_id)
raw.tcga.count<-assays(rse)[[1]]
colnames(raw.tcga.count)<-rse@colData@listData$gdc_cases.submitter_id
length(rse@rowRanges$symbol)
index<- match(toupper(genes),toupper(rse@rowRanges$symbol@listData))
test_index<-na.omit(index)
length(test_index)
sum(is.null(rse@rowRanges$symbol@listData[test_index]))
raw.tcga.count.reduced<-raw.tcga.count[as.numeric(test_index),]
rownames(raw.tcga.count.reduced)<-rse@rowRanges$symbol@listData[test_index]
dim(raw.tcga.count.reduced)
tcga_genes<-rse@rowRanges$symbol@listData[test_index]
rse_reduced <- rse@colData@listData$gdc_cases.samples.sample_type == 'Primary Tumor' & rse@colData@listData$gdc_cases.project.project_id == 'TCGA-GBM'
sum(rse_reduced)

# 157 Primary GBM Tumor Data
raw.tcga.count.reduced.primary <- raw.tcga.count.reduced[,rse_reduced]
rownames(raw.tcga.count.reduced.primary)
colnames(raw.tcga.count.reduced.primary)
dim(na.omit(raw.tcga.count.reduced.primary))
# Dimension 1888*157
save(raw.tcga.count.reduced.primary,file="gbm_primary_tumor_expression_data.Rdata")

# Clincial Data of age,gender,class of tumor,Karnofsky Score
age_at_diagnosis=rse@colData@listData$gdc_cases.diagnoses.age_at_diagnosis/365
length(age_at_diagnosis)
age<-age_at_diagnosis[rse_reduced]
sum(is.na(age))
length(age)
gender_all<-rse@colData@listData$gdc_cases.demographic.gender
gender<-gender_all[rse_reduced]
sum(is.na(gender))
length(gender)
score_all<-rse$xml_karnofsky_performance_score
score<-score_all[rse_reduced]
sum(is.na(score))
length(score)
clinical_data<- data.frame('feature'=colnames(raw.tcga.count.reduced.primary),'age'=age,'gender'=gender,'score'=score)
save(clinical_data,file="gbm_clinical_data.Rdata")

# Creating survival data  
follow_up_days=rse@colData@listData$gdc_cases.diagnoses.days_to_last_follow_up
length(follow_up_days)
follow_up_days.primary <- follow_up_days[rse_reduced]
length(follow_up_days.primary)
# Creating survival data  
follow_up_days=rse@colData@listData$gdc_cases.diagnoses.days_to_last_follow_up
length(follow_up_days)
follow_up_days.primary <- follow_up_days[rse_reduced]
length(follow_up_days.primary)

days_to_death=rse@colData@listData$gdc_cases.diagnoses.days_to_death
length(days_to_death)
days_to_death.primary <- days_to_death[rse_reduced]
length(days_to_death.primary)


die=rse@colData@listData$gdc_cases.diagnoses.vital_status
length(die)
die.primary <- die[rse_reduced]
length(die.primary)

event_days<-list()
for (i in c(1:length(follow_up_days.primary))) {
  if (is.na(follow_up_days.primary[i])){
    #print(is.na(days_to_death.primary[i]))
    event_days[i]=days_to_death.primary[i]
  }
  else{
    event_days[i]=follow_up_days.primary[i]
  }
}
survival_data<- data.frame('feature'=colnames(raw.tcga.count.reduced.primary),'time'=unlist(event_days),'status'=die.primary)
save(survival_data,file="gbm_survival_data.Rdata")

# Keeping Data which has all the information
# load the expression,clinical,survival data
rm(list=ls(all=TRUE))
load('gbm_primary_tumor_expression_data.Rdata')
load('gbm_clinical_data.Rdata')
load('gbm_survival_data.Rdata')
dim(na.omit(raw.tcga.count.reduced.primary))
dim(unique(na.omit(clinical_data)))
dim(na.omit(survival_data))
# Removing the NA from everywhere
col_id<-unique(na.omit(clinical_data)$feature)
survival_data.reduced <- survival_data[match(col_id,survival_data$feature),]
final_col_id<-na.omit(survival_data.reduced)$feature
length(final_col_id)
# Final Dataset
gbm_mRNA <- raw.tcga.count.reduced.primary[,final_col_id]
dim(gbm_mRNA)
gbm_clincial<- clinical_data[match(final_col_id,clinical_data$feature),]
dim(gbm_clincial)
gbm_survival <- survival_data[match(final_col_id,survival_data$feature),]
dim(gbm_survival)
col_id<-unique(na.omit(clinical_data)$feature)
survival_data.reduced <- survival_data[match(col_id,survival_data$feature),]
final_col_id<-na.omit(survival_data.reduced)$feature
length(final_col_id)
# Final Dataset
gbm_mRNA <- raw.tcga.count.reduced.primary[,final_col_id]
dim(gbm_mRNA)
gbm_clincial<- clinical_data[match(final_col_id,clinical_data$feature),]
dim(gbm_clincial)
gbm_survival <- survival_data[match(final_col_id,survival_data$feature),]
dim(gbm_survival)

save(gbm_mRNA,file = 'gbm_final_mRNA_raw_count.data')
save(gbm_clincial,file = 'gbm_final_clincial.data')
save(gbm_survival,file = 'gbm_final_survival.data')

#---
# LUSC Data preprocessing
#---
#Clear workspace
rm(list=ls(all=TRUE))
load('only_ctc_updated_filter.Rdata')
dim(filter_data)
#load Lung dataset 
load('rse_gene_lung.Rdata')
dim(rse_gene)
rse<-scale_counts(rse_gene)
rm(rse_gene)
dim(rse)
table(rse@colData@listData$gdc_cases.project.project_id)
raw.tcga.count<-assays(rse)[[1]]
colnames(raw.tcga.count)<-rse@colData@listData$gdc_cases.submitter_id
length(rse@rowRanges$symbol)
index<- match(toupper(rownames(filter_data)),rse@rowRanges$symbol@listData)
test_index<-na.omit(index)
length(test_index)
raw.tcga.count.reduced<-raw.tcga.count[as.numeric(test_index),]
rownames(raw.tcga.count.reduced)<-rse@rowRanges$symbol@listData[test_index]
rownames(raw.tcga.count.reduced)
colnames(raw.tcga.count.reduced)
dim(raw.tcga.count.reduced)
#dim(rse)
rse_reduced <- rse@colData@listData$gdc_cases.samples.sample_type == 'Primary Tumor' & rse@colData@listData$gdc_cases.project.project_id == 'TCGA-LUSC'
sum(rse_reduced)

# 504 Primary Lung Tumor Data
raw.tcga.count.reduced.primary <- raw.tcga.count.reduced[,rse_reduced]
rownames(raw.tcga.count.reduced.primary)
colnames(raw.tcga.count.reduced.primary)
dim(na.omit(raw.tcga.count.reduced.primary))
# Dimension 1873*504
save(raw.tcga.count.reduced.primary,file="lung_primary_tumor_expression_data.Rdata")

# Clincial Data of age,gender,stage of tumor,Karnofsky Score
age_at_diagnosis=rse@colData@listData$gdc_cases.diagnoses.age_at_diagnosis/365
length(age_at_diagnosis)
age<-age_at_diagnosis[rse_reduced]
sum(is.na(age))
length(age)
gender_all<-rse@colData@listData$gdc_cases.demographic.gender
gender<-gender_all[rse_reduced]
sum(is.na(gender))
length(gender)
# score_all<-rse$xml_karnofsky_performance_score
# score<-score_all[rse_reduced]
# sum(is.na(score))
# length(score)
cancer_stage_all<-rse$xml_stage_event_pathologic_stage
cancer_stage<-cancer_stage_all[rse_reduced]
sum(is.na(cancer_stage))
length(cancer_stage)
clinical_data<- data.frame('feature'=colnames(raw.tcga.count.reduced.primary),'age'=age,'gender'=gender,'stage'=cancer_stage)
save(clinical_data,file="lung_clinical_data.Rdata")

# Creating survival data  
follow_up_days=rse@colData@listData$gdc_cases.diagnoses.days_to_last_follow_up
length(follow_up_days)
follow_up_days.primary <- follow_up_days[rse_reduced]
length(follow_up_days.primary)

days_to_death=rse@colData@listData$gdc_cases.diagnoses.days_to_death
length(days_to_death)
days_to_death.primary <- days_to_death[rse_reduced]
length(days_to_death.primary)


die=rse@colData@listData$gdc_cases.diagnoses.vital_status
length(die)
die.primary <- die[rse_reduced]
length(die.primary)

event_days<-list()
for (i in c(1:length(follow_up_days.primary))) {
  if (is.na(follow_up_days.primary[i])){
    #print(is.na(days_to_death.primary[i]))
    event_days[i]=days_to_death.primary[i]
  }
  else{
    event_days[i]=follow_up_days.primary[i]
  }
}
survival_data<- data.frame('feature'=colnames(raw.tcga.count.reduced.primary),'time'=unlist(event_days),'status'=die.primary)
save(survival_data,file="lung_survival_data.Rdata")

#Clear workspace
#rm(list=ls(all=TRUE))
# Keeping Data which has all the information
# load the expression,clinical,survival data
load('lung_primary_tumor_expression_data.Rdata')
load('lung_clinical_data.Rdata')
load('lung_survival_data.Rdata')
dim(na.omit(raw.tcga.count.reduced.primary))
dim(na.omit(clinical_data))
dim(na.omit(survival_data))
# Removing the NA from everywhere
col_id<-na.omit(clinical_data)$feature
survival_data.reduced <- survival_data[match(col_id,survival_data$feature),]
final_col_id<-na.omit(survival_data.reduced)$feature
length(final_col_id)
length(unique(final_col_id))
# Final Dataset
lung_mRNA <- raw.tcga.count.reduced.primary[,unique(final_col_id)]
dim(lung_mRNA)
lung_clincial<- clinical_data[match(unique(final_col_id),clinical_data$feature),]
dim(lung_clincial)
lung_survival <- survival_data[match(unique(final_col_id),survival_data$feature),]
dim(lung_survival)
save(lung_mRNA,file = 'lung_final_mRNA_raw_count.data')
save(lung_clincial,file = 'lung_final_clincial.data')
save(lung_survival,file = 'lung_final_survival.data')


#Clear workspace
#rm(list=ls(all=TRUE))
load('only_ctc_updated_filter.Rdata')
dim(filter_data)
#load OV dataset 
load('rse_gene_ovary.Rdata')
dim(rse_gene)
rse<-scale_counts(rse_gene)
rm(rse_gene)
dim(rse)
table(rse@colData@listData$gdc_cases.project.project_id)
raw.tcga.count<-assays(rse)[[1]]
colnames(raw.tcga.count)<-rse@colData@listData$gdc_cases.submitter_id
length(rse@rowRanges$symbol)
index<- match(toupper(rownames(filter_data)),rse@rowRanges$symbol@listData)
test_index<-na.omit(index)
length(test_index)
raw.tcga.count.reduced<-raw.tcga.count[as.numeric(test_index),]
rownames(raw.tcga.count.reduced)<-rse@rowRanges$symbol@listData[test_index]
rownames(raw.tcga.count.reduced)
colnames(raw.tcga.count.reduced)
dim(raw.tcga.count.reduced)
#dim(rse)
rse_reduced <- rse@colData@listData$gdc_cases.samples.sample_type == 'Primary Tumor'
sum(rse_reduced)

# 422 Primary OV Tumor Data
raw.tcga.count.reduced.primary <- raw.tcga.count.reduced[,rse_reduced]
rownames(raw.tcga.count.reduced.primary)
colnames(raw.tcga.count.reduced.primary)
dim(na.omit(raw.tcga.count.reduced.primary))
# Dimension 1873*422
save(raw.tcga.count.reduced.primary,file="ov_primary_tumor_expression_data.Rdata")

# Clincial Data of age,gender,stage of tumo/clincial_stage,Karnofsky Score
age_at_diagnosis=rse@colData@listData$gdc_cases.diagnoses.age_at_diagnosis/365
length(age_at_diagnosis)
age<-age_at_diagnosis[rse_reduced]
sum(is.na(age))
length(age)
gender_all<-rse@colData@listData$gdc_cases.demographic.gender
gender<-gender_all[rse_reduced]
sum(is.na(gender))
length(gender)
# score_all<-rse$xml_karnofsky_performance_score
# score<-score_all[rse_reduced]
# sum(is.na(score))
# length(score)
cancer_stage_all<-rse$xml_stage_event_clinical_stage
cancer_stage<-cancer_stage_all[rse_reduced]
sum(is.na(cancer_stage))
length(cancer_stage)
clinical_data<- data.frame('feature'=colnames(raw.tcga.count.reduced.primary),'gender'=gender,'age'=age,'stage'=cancer_stage)
save(clinical_data,file="ov_clinical_data.Rdata")

# Creating survival data  
follow_up_days=rse@colData@listData$gdc_cases.diagnoses.days_to_last_follow_up
length(follow_up_days)
follow_up_days.primary <- follow_up_days[rse_reduced]
length(follow_up_days.primary)

days_to_death=rse@colData@listData$gdc_cases.diagnoses.days_to_death
length(days_to_death)
days_to_death.primary <- days_to_death[rse_reduced]
length(days_to_death.primary)


die=rse@colData@listData$gdc_cases.diagnoses.vital_status
length(die)
die.primary <- die[rse_reduced]
length(die.primary)

event_days<-list()
for (i in c(1:length(follow_up_days.primary))) {
  if (is.na(follow_up_days.primary[i])){
    #print(is.na(days_to_death.primary[i]))
    event_days[i]=days_to_death.primary[i]
  }
  else{
    event_days[i]=follow_up_days.primary[i]
  }
}
survival_data<- data.frame('feature'=colnames(raw.tcga.count.reduced.primary),'time'=unlist(event_days),'status'=die.primary)
save(survival_data,file="ov_survival_data.Rdata")




#Clear workspace
rm(list=ls(all=TRUE))
# Keeping Data which has all the information
# load the expression,clinical,survival data
load('ov_primary_tumor_expression_data.Rdata')
load('ov_clinical_data.Rdata')
load('ov_survival_data.Rdata')
dim(na.omit(raw.tcga.count.reduced.primary))
dim(na.omit(clinical_data))
dim(na.omit(survival_data))
# Removing the NA from everywhere
col_id<-na.omit(clinical_data)$feature
survival_data.reduced <- survival_data[match(col_id,survival_data$feature),]
final_col_id<-na.omit(survival_data.reduced)$feature
length(final_col_id)
length(unique(final_col_id))
# Final Dataset
ov_mRNA <- raw.tcga.count.reduced.primary[,unique(final_col_id)]
dim(ov_mRNA)
ov_clincial<- clinical_data[match(unique(final_col_id),clinical_data$feature),]
dim(ov_clincial)
ov_survival <- survival_data[match(unique(final_col_id),survival_data$feature),]
dim(ov_survival)
save(ov_mRNA,file = 'ov_final_mRNA_raw_count.data')
save(ov_clincial,file = 'ov_final_clincial.data')
save(ov_survival,file = 'ov_final_survival.data')


#Clear workspace
rm(list=ls(all=TRUE))
load('only_ctc_updated_filter.Rdata')
dim(filter_data)
#load KIRC dataset 
load('rse_gene_kidney.Rdata')
dim(rse_gene)
rse<-scale_counts(rse_gene)
rm(rse_gene)
dim(rse)
table(rse@colData@listData$gdc_cases.project.project_id)
raw.tcga.count<-assays(rse)[[1]]
colnames(raw.tcga.count)<-rse@colData@listData$gdc_cases.submitter_id
length(rse@rowRanges$symbol)
index<- match(toupper(rownames(filter_data)),rse@rowRanges$symbol@listData)
test_index<-na.omit(index)
length(test_index)
raw.tcga.count.reduced<-raw.tcga.count[as.numeric(test_index),]
rownames(raw.tcga.count.reduced)<-rse@rowRanges$symbol@listData[test_index]
rownames(raw.tcga.count.reduced)
colnames(raw.tcga.count.reduced)
dim(raw.tcga.count.reduced)
#dim(rse)
rse_reduced <- rse@colData@listData$gdc_cases.samples.sample_type == 'Primary Tumor' & rse@colData@listData$gdc_cases.project.project_id == 'TCGA-KIRC'
sum(rse_reduced)

# 543 Primary Kidney Tumor Data
raw.tcga.count.reduced.primary <- raw.tcga.count.reduced[,rse_reduced]
rownames(raw.tcga.count.reduced.primary)
colnames(raw.tcga.count.reduced.primary)
dim(na.omit(raw.tcga.count.reduced.primary))
# Dimension 1873*543
save(raw.tcga.count.reduced.primary,file="kidney_primary_tumor_expression_data.Rdata")

# Clincial Data of age,gender,stage of tumo/clincial_stage,Karnofsky Score
age_at_diagnosis=rse@colData@listData$gdc_cases.diagnoses.age_at_diagnosis/365
length(age_at_diagnosis)
age<-age_at_diagnosis[rse_reduced]
sum(is.na(age))
length(age)
gender_all<-rse@colData@listData$gdc_cases.demographic.gender
gender<-gender_all[rse_reduced]
sum(is.na(gender))
length(gender)
# score_all<-rse$xml_karnofsky_performance_score
# score<-score_all[rse_reduced]
# sum(is.na(score))
# length(score)
cancer_stage_all<-rse$xml_stage_event_pathologic_stage
cancer_stage<-cancer_stage_all[rse_reduced]
sum(is.na(cancer_stage))
length(cancer_stage)

clinical_data<- data.frame('feature'=colnames(raw.tcga.count.reduced.primary),'gender'=gender,'age'=age,'stage'=cancer_stage)
save(clinical_data,file="kidney_clinical_data.Rdata")

# Creating survival data  
follow_up_days=rse@colData@listData$gdc_cases.diagnoses.days_to_last_follow_up
length(follow_up_days)
follow_up_days.primary <- follow_up_days[rse_reduced]
length(follow_up_days.primary)

days_to_death=rse@colData@listData$gdc_cases.diagnoses.days_to_death
length(days_to_death)
days_to_death.primary <- days_to_death[rse_reduced]
length(days_to_death.primary)


die=rse@colData@listData$gdc_cases.diagnoses.vital_status
length(die)
die.primary <- die[rse_reduced]
length(die.primary)

event_days<-list()
for (i in c(1:length(follow_up_days.primary))) {
  if (is.na(follow_up_days.primary[i])){
    #print(is.na(days_to_death.primary[i]))
    event_days[i]=days_to_death.primary[i]
  }
  else{
    event_days[i]=follow_up_days.primary[i]
  }
}
survival_data<- data.frame('feature'=colnames(raw.tcga.count.reduced.primary),'time'=unlist(event_days),'status'=die.primary)
save(survival_data,file="kidney_survival_data.Rdata")


#Clear workspace
rm(list=ls(all=TRUE))
# Keeping Data which has all the information
# load the expression,clinical,survival data
load('kidney_primary_tumor_expression_data.Rdata')
load('kidney_clinical_data.Rdata')
load('kidney_survival_data.Rdata')
dim(na.omit(raw.tcga.count.reduced.primary))
dim(na.omit(clinical_data))
dim(na.omit(survival_data))
# Removing the NA from everywhere
col_id<-na.omit(clinical_data)$feature
survival_data.reduced <- survival_data[match(col_id,survival_data$feature),]
final_col_id<-na.omit(survival_data.reduced)$feature
length(final_col_id)
length(unique(final_col_id))
# Final Dataset
kidney_mRNA <- raw.tcga.count.reduced.primary[,unique(final_col_id)]
dim(kidney_mRNA)
kidney_clincial<- clinical_data[match(unique(final_col_id),clinical_data$feature),]
dim(kidney_clincial)
kidney_survival <- survival_data[match(unique(final_col_id),survival_data$feature),]
dim(kidney_survival)
save(kidney_mRNA,file = 'kidney_final_mRNA_raw_count.data')
save(kidney_clincial,file = 'kidney_final_clincial.data')
save(kidney_survival,file = 'kidney_final_survival.data')