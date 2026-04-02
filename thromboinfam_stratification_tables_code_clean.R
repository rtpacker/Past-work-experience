# =============================================================================
# Residential Segregation & Thrombo-Inflammatory Biomarkers
# Stratification Tables — Survey-Weighted Models (svyglm)
# =============================================================================
# Description: Fits survey-weighted linear models within race/gender strata
#              and exports results to Excel workbooks (one per biomarker
#              per stratum).
#
# Strata:      White, Black, Male, Female,
#              Black Male, Black Female, White Male, White Female
# Exposures:   Dissimilarity index, Interaction index, Isolation index (std)
# Outcomes:    log_ddimer, log_crp, log_ifn, log_tnf, log_il6, eselectin, fix
# Models:      Unadjusted, Age-adjusted
#
# 
# =============================================================================

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

library(knitr)
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




# Create new standardized variables  "rs_dissimilarityb", "rs_interactionb", "rs_isolationb"
analyticsamplethr$dissimilarity_index_std <- analyticsamplethr$rs_dissimilarityb/ sd(analyticsamplethr$rs_dissimilarityb)
analyticsamplethr$interaction_index_std <- analyticsamplethr$rs_interactionb / sd(analyticsamplethr$rs_interactionb)
analyticsamplethr$isolation_index_std <- analyticsamplethr$rs_isolationb / sd(analyticsamplethr$rs_isolationb)


# #make tertitles
# analyticsamplethr<-analyticsamplethr %>%mutate(dissimilarity_index_tertiles = ntile(rs_dissimilarityb, 3))
# analyticsamplethr<-analyticsamplethr %>%mutate(interaction_index_tertiles = ntile(rs_interactionb, 3))
# analyticsamplethr<-analyticsamplethr %>%mutate(isolation_index_tertiles = ntile(rs_isolationb, 3))



# #TERTILES FOR THE INDINCES
# analyticsamplethr$dissimilarity_index_tertiles <- factor(analyticsamplethr$dissimilarity_index_tertiles, levels = c("1", "2","3"))
# analyticsamplethr$dissimilarity_index_tertiles <- relevel(analyticsamplethr$dissimilarity_index_tertiles, ref = "1")
# analyticsamplethr$interaction_index_tertiles <- factor(analyticsamplethr$interaction_index_tertiles, levels = c("1", "2","3"))
# analyticsamplethr$interaction_index_tertiles <- relevel(analyticsamplethr$interaction_index_tertiles, ref = "1")
# analyticsamplethr$isolation_index_tertiles <- factor(analyticsamplethr$isolation_index_tertiles, levels = c("1", "2","3"))
# analyticsamplethr$isolation_index_tertiles <- relevel(analyticsamplethr$isolation_index_tertiles, ref = "1")





#Covariates:age, gender, race, income_4cat, ed_cat
#nses_tract, alc_niaaa,smoke, exercise_cat,diet7(LS7dietscore)
##diab_srmed_glu, stroke_sr,tia_sr, cad_sr_ecg 



#########GENDER SUBSETS
analyticsamplethr_GSTRAT_M <- analyticsamplethr %>% filter(gender == "M")
analyticsamplethr_GSTRAT_F=analyticsamplethr%>%filter(gender=="F")

analyticsamplethr_RSTRAT_W=analyticsamplethr%>%filter(race=="W")
analyticsamplethr_RSTRAT_B=analyticsamplethr%>%filter(race=="B")



analyticsamplethr_GRSTRAT_MW=analyticsamplethr%>%filter(gender=="M",race=="W")
analyticsamplethr_GRSTRAT_MB=analyticsamplethr%>%filter(gender=="M",race=="B")

analyticsamplethr_GRSTRAT_FW=analyticsamplethr%>%filter(gender=="F",race=="W")
analyticsamplethr_GRSTRAT_FB=analyticsamplethr%>%filter(gender=="F",race=="B")




#Covariates:age, gender, race, income_4cat, ed_cat
#nses_tract, alc_niaaa,smoke, exercise_cat,diet7(LS7dietscore)
##diab_srmed_glu, stroke_sr,tia_sr, cad_sr_ecg 

w_dgn<- svydesign(~1, prob = ~prob,
                           strata=~race*gender,
                           data = analyticsamplethr %>%
                             filter(race=="W")) 

b_dgn<- svydesign(~1, prob = ~prob,
                  strata=~race*gender,
                  data = analyticsamplethr %>%
                    filter(race=="B")) 

men_dgn <- svydesign(~1, prob = ~prob,
                     strata=~race*gender,
                     data = analyticsamplethr %>%
                       filter(gender=="M"))                                 

women_dgn <- svydesign(~1, prob = ~prob,
                       strata=~race*gender,
                       data = analyticsamplethr %>%
                         filter(gender=="F"))


black_male_dgn <- svydesign(~1, prob = ~prob, strata=~race*gender, data = analyticsamplethr %>% filter(race=="B" & gender=="M"))
black_female_dgn <- svydesign(~1, prob = ~prob, strata=~race*gender, data = analyticsamplethr %>% filter(race=="B" & gender=="F"))
white_male_dgn <- svydesign(~1, prob = ~prob, strata=~race*gender, data = analyticsamplethr %>% filter(race=="W" & gender=="M"))
white_female_dgn <- svydesign(~1, prob = ~prob, strata=~race*gender, data = analyticsamplethr %>% filter(race=="W" & gender=="F"))


#####

#ddimer isolation and interaction gender
#CRP isolation and interaction gender and race
#eselectin

# List of exposures
exposures <- c("dissimilarity_index_std", "interaction_index_std", "isolation_index_std")  # exposures are the 3 index cats
b <- c("log_ddimer", "eselectin","log_ifn","log_crp")#makesure fix isnt log
# List of models
models <- c("Unadjusted", "Age and gender adjusted")
# List to store model results
model_results <- list()

# Iterate over exposures
for (exposure in exposures) {
  # Iterate over biomarkers
  for (biomarker in b) {
    # Formulate the formula for each exposure and biomarker
    formula_unadjusted <- as.formula(paste(biomarker, "~", exposure))
    formula_age_gender <- as.formula(paste(biomarker, "~", exposure, "+ age"))
    
    
    
    
    
    
    # Fit the models
    model_unadjusted <- svyglm(formula_unadjusted,design=black_male_dgn, data = analyticsamplethr_GRSTRAT_MB)
    model_age_gender <- svyglm(formula_age_gender,design=black_male_dgn, data = analyticsamplethr_GRSTRAT_MB)
    
    
    # Store the model results
    model_results[[paste(exposure, "_", biomarker, sep = "")]] <- list(
      Unadjusted = summary(model_unadjusted),
      `Age and gender adjusted` = summary(model_age_gender)
      
      
    )
  }
}

#LOOP WORKS
#Need to do: standardize...Fix model results output table(kat c code)....log the biomarkers

# Calculate weighted standard deviation for blineFVIII
# loweprojanset_exlu_t$blineFVIII_sd <-  wtd.var(loweprojanset_exlu_t$blineFVIII, weights = loweprojanset_exlu_t$cc_weight)
# loweprojanset_exlu_t$blineFVIII_sd <- sqrt(loweprojanset_exlu_t$blineFVIII_sd)

#create 3 new rs vars that is the index/sd(index) and use those ij the models instead!

###REMOVE THE REVERSE T AND REVERSET CI for fix and eselectin add in normal CI 95% 


#fix######reverse transform for this is wonky cuz this bio makrer is not logged

models_of_interest <- c("Unadjusted", "Age and gender adjusted")

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$dissimilarity_index_std_log_crp[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err <- model$coefficients["dissimilarity_index_std", "Std. Error"]
  p_value <- model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df dataframe
  results_df <- rbind(results_df, result)
}

# Write dissimilarity results to the second sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)

# Print or further process results_df as needed
print(results_df)



for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$isolation_index_std_log_crp[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["isolation_index_std", "Estimate"]
  std_err <- model$coefficients["isolation_index_std", "Std. Error"]
  p_value <- model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  if (lcl_p > ucl_p) {
    # Swap values if lcl_p is greater than ucl_p
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 dataframe
  results_df_2 <- rbind(results_df_2, result)
}

# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2, colNames = TRUE)



# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$interaction_index_std_log_crp[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["interaction_index_std", "Estimate"]
  std_err <- model$coefficients["interaction_index_std", "Std. Error"]
  p_value <- model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 dataframe
  results_df_3 <- rbind(results_df_3, result)
}

# Write interaction results to the second sheet
writeData(wb, sheet = "Interaction_Results", results_df_3, colNames = TRUE)

# Print or further process results_df_3 as needed
print(results_df_3)



# Save the Excel workbook
saveWorkbook(wb, "path/to/output/crpRaceBGenMstrat2.xlsx")
























