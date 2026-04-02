# =============================================================================
# Residential Segregation & Thrombo-Inflammatory Biomarkers
# Interaction Tables — Survey-Weighted Models (svyglm)
# =============================================================================
# Description: Tests for interaction between residential segregation indices
#              and race/gender on thrombo-inflammatory biomarkers using
#              survey-weighted linear models (svyglm).
#
# Exposures:   Dissimilarity index, Interaction index, Isolation index (std)
# Outcomes:    log_ddimer, log_crp, log_ifn, log_tnf, log_il6, eselectin, fix
# Interactions tested: x Race, x Gender, x Race*Gender
#

# =============================================================================




###############

###############
library(dplyr)
library(data.table)
library(zoo) 
library(tidyr)
library(lubridate) 
library(tidyverse)
library(knitr)
library(kableExtra)
library(sas7bdat)
library(haven)
library("nimble")
library(xlsx) 
library(openxlsx)
library(survey)
library(broom)
library(purrr)


#read in data..... IL-1b needs to be broken into quartiles (see analytic plan)

# #il1b in this data.....il1....confirm but combed through other data and this look like it....need to be in quartiles
# biomedior=read.sas7bdat('path/to/biomedior_results_id.sas7bdat')
# 
# 
# BMDR_dat=readRDS('path/to/BMDR_dat.rds')
# calvars_rsi=read.sas7bdat('path/to/calcvars_rsi.sas7bdat')


#analytic sample data


analyticsamplethr=readRDS('path/to/analyticsamplethr.rds')
head(analyticsamplethr)#4362 before exclusions....2534 after exclusions


analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr$rs_dissimilarityb), ]#none
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr$rs_interactionb), ]#none
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr$rs_isolationb), ]#none

analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_ddimer"]]), ]#148
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_crp"]]), ]#48
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_ifn"]]), ]#37
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_tnf"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_il6"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["eselectin"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["fix"]]), ]#6
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["il1"]]), ]#6


analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["age"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["gender"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["race"]]), ]#0


analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_ddimer"]]), ]#148
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_crp"]]), ]#48
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_ifn"]]), ]#37
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_tnf"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_il6"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["eselectin"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["fix"]]), ]#6
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["il1"]]), ]#6

analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["age"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["gender"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["race"]]), ]#0

# Create new standardized variables  "rs_dissimilarityb", "rs_interactionb", "rs_isolationb"
analyticsamplethr$dissimilarity_index_std <- analyticsamplethr$rs_dissimilarityb/ sd(analyticsamplethr$rs_dissimilarityb)
analyticsamplethr$interaction_index_std <- analyticsamplethr$rs_interactionb / sd(analyticsamplethr$rs_interactionb)
analyticsamplethr$isolation_index_std <- analyticsamplethr$rs_isolationb / sd(analyticsamplethr$rs_isolationb)


#make tertitles
analyticsamplethr<-analyticsamplethr %>%mutate(dissimilarity_index_tertiles = ntile(rs_dissimilarityb, 3))
analyticsamplethr<-analyticsamplethr %>%mutate(interaction_index_tertiles = ntile(rs_interactionb, 3))
analyticsamplethr<-analyticsamplethr %>%mutate(isolation_index_tertiles = ntile(rs_isolationb, 3))




#TERTILES FOR THE INDINCES
analyticsamplethr$dissimilarity_index_tertiles <- factor(analyticsamplethr$dissimilarity_index_tertiles, levels = c("1", "2","3"))
analyticsamplethr$dissimilarity_index_tertiles <- relevel(analyticsamplethr$dissimilarity_index_tertiles, ref = "1")
analyticsamplethr$interaction_index_tertiles <- factor(analyticsamplethr$interaction_index_tertiles, levels = c("1", "2","3"))
analyticsamplethr$interaction_index_tertiles <- relevel(analyticsamplethr$interaction_index_tertiles, ref = "1")
analyticsamplethr$isolation_index_tertiles <- factor(analyticsamplethr$isolation_index_tertiles, levels = c("1", "2","3"))
analyticsamplethr$isolation_index_tertiles <- relevel(analyticsamplethr$isolation_index_tertiles, ref = "1")







dgn= svydesign(~1, prob = ~prob, strata=~race*gender,data = analyticsamplethr)




hist(analyticsamplethr$eselectin)


exposures <- c("dissimilarity_index_std", "interaction_index_std", "isolation_index_std","isolation_index_tertiles","dissimilarity_index_tertiles","interaction_index_tertiles")  # exposures are the 3 index cats
b <- c("log_ddimer", "log_crp", "log_ifn", "log_tnf", "log_il6", "eselectin","fix")#makesure fix isnt log







#######NOTE: to get this code to work...just copy and replace the exposures and cycle through below code...copy and pate table results into excel





#fix
# Fit the unadjusted model
model_unadjusted <- svyglm(fix ~ interaction_index_std + age + gender,design=dgn, data = analyticsamplethr)


# Assess the interaction with race
interaction_race <- svyglm(fix ~ interaction_index_std * race + age + gender + race,design=dgn, data = analyticsamplethr)
summary(interaction_race)

# Assess the interaction with gender
interaction_gender <- svyglm(fix ~ interaction_index_std * gender + age + race +gender,design=dgn, data = analyticsamplethr)
summary(interaction_gender)

# Assess the interaction with both race and gender
interaction_race_gender <- svyglm(fix ~ interaction_index_std * race * gender + age +race+gender,design=dgn, data = analyticsamplethr)
summary(interaction_race_gender)




#ddimer
# Fit the unadjusted model
model_unadjusted <- svyglm(log_ddimer ~ interaction_index_std + age + gender,design=dgn, data = analyticsamplethr)


# Assess the interaction with race
interaction_race <- svyglm(log_ddimer ~ interaction_index_std * race + age + gender + race,design=dgn, data = analyticsamplethr)
summary(interaction_race)

# Assess the interaction with gender
interaction_gender <- svyglm(log_ddimer ~ interaction_index_std * gender + age + race +gender,design=dgn, data = analyticsamplethr)
summary(interaction_gender)

# Assess the interaction with both race and gender
interaction_race_gender <- svyglm(log_ddimer ~ interaction_index_std * race * gender + age +race+gender,design=dgn, data = analyticsamplethr)
summary(interaction_race_gender)


#crp
# Fit the unadjusted model
model_unadjusted <- svyglm(log_crp ~ interaction_index_std + age + gender,design=dgn, data = analyticsamplethr)


# Assess the interaction with race
interaction_race <- svyglm(log_crp ~ interaction_index_std * race + age + gender + race,design=dgn, data = analyticsamplethr)
summary(interaction_race)

# Assess the interaction with gender
interaction_gender <- svyglm(log_crp ~ interaction_index_std * gender + age + race +gender,design=dgn, data = analyticsamplethr)
summary(interaction_gender)

# Assess the interaction with both race and gender
interaction_race_gender <- svyglm(log_crp ~ interaction_index_std * race * gender + age +race+gender,design=dgn, data = analyticsamplethr)
summary(interaction_race_gender)


#ifn
# Fit the unadjusted model
model_unadjusted <- svyglm(log_ifn ~ interaction_index_std + age + gender,design=dgn, data = analyticsamplethr)


# Assess the interaction with race
interaction_race <- svyglm(log_ifn ~ interaction_index_std * race + age + gender + race,design=dgn, data = analyticsamplethr)
summary(interaction_race)

# Assess the interaction with gender
interaction_gender <- svyglm(log_ifn ~ interaction_index_std * gender + age + race +gender,design=dgn, data = analyticsamplethr)
summary(interaction_gender)

# Assess the interaction with both race and gender
interaction_race_gender <- svyglm(log_ifn ~ interaction_index_std * race * gender + age +race+gender,design=dgn, data = analyticsamplethr)
summary(interaction_race_gender)



#tnf
# Fit the unadjusted model
model_unadjusted <- svyglm(log_tnf ~ interaction_index_std + age + gender,design=dgn, data = analyticsamplethr)


# Assess the interaction with race
interaction_race <- svyglm(log_tnf ~ interaction_index_std * race + age + gender + race,design=dgn, data = analyticsamplethr)
summary(interaction_race)

# Assess the interaction with gender
interaction_gender <- svyglm(log_tnf ~ interaction_index_std * gender + age + race +gender,design=dgn, data = analyticsamplethr)
summary(interaction_gender)

# Assess the interaction with both race and gender
interaction_race_gender <- svyglm(log_tnf ~ interaction_index_std * race * gender + age +race+gender,design=dgn, data = analyticsamplethr)
summary(interaction_race_gender)


#il6
# Fit the unadjusted model
model_unadjusted <- svyglm(log_il6 ~ interaction_index_std + age + gender,design=dgn, data = analyticsamplethr)


# Assess the interaction with race
interaction_race <- svyglm(log_il6 ~ interaction_index_std * race + age + gender + race,design=dgn, data = analyticsamplethr)
summary(interaction_race)

# Assess the interaction with gender
interaction_gender <- svyglm(log_il6 ~ interaction_index_std * gender + age + race +gender,design=dgn, data = analyticsamplethr)
summary(interaction_gender)

# Assess the interaction with both race and gender
interaction_race_gender <- svyglm(log_il6 ~ interaction_index_std * race * gender + age +race+gender,design=dgn, data = analyticsamplethr)
summary(interaction_race_gender)




#eselectin
# Fit the unadjusted model
model_unadjusted <- svyglm(eselectin ~ interaction_index_std + age + gender,design=dgn, data = analyticsamplethr)


# Assess the interaction with race
interaction_race <- svyglm(eselectin ~ interaction_index_std * race + age + gender + race,design=dgn, data = analyticsamplethr)
summary(interaction_race)

# Assess the interaction with gender
interaction_gender <- svyglm(eselectin ~ interaction_index_std * gender + age + race +gender,design=dgn, data = analyticsamplethr)
summary(interaction_gender)

# Assess the interaction with both race and gender
interaction_race_gender <- svyglm(eselectin ~ interaction_index_std * race * gender + age +race+gender,design=dgn, data = analyticsamplethr)
summary(interaction_race_gender)

# #Il1
# #isolation_index_tertiles_log_il1,dissimilarity_index_tertiles_log_il1,interaction_index_tertiles_log_il1
# # Fit the unadjusted model
# model_unadjusted <- svyglm(log_il1 ~ dissimilarity_index_tertiles + age + gender,design=dgn, data = analyticsamplethr)
# 
# 
# # Assess the interaction with race
# interaction_race <- svyglm(log_il1 ~ dissimilarity_index_tertiles * race + age + gender + race,design=dgn, data = analyticsamplethr)
# summary(interaction_race)
# 
# # Assess the interaction with gender
# interaction_gender <- svyglm(log_il1 ~ dissimilarity_index_tertiles * gender + age + race +gender,design=dgn, data = analyticsamplethr)
# summary(interaction_gender)
# 
# # Assess the interaction with both race and gender
# interaction_race_gender <- svyglm(log_il1 ~ dissimilarity_index_tertiles * race * gender + age +race+gender,design=dgn, data = analyticsamplethr)
# summary(interaction_race_gender)
