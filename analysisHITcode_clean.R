# =============================================================================
# HIT (Heparin-Induced Thrombocytopenia) Analysis
# =============================================================================
# Description: Baseline characteristics table and univariate logistic regression
#              models examining factors associated with anticoagulation agent
#              selection following a positive HIT antibody lab result.
#
# Outcome:     outcome_dayafterhitlabacagent (binary: 1/0)
# Data:        Multi-site analytic dataset (hitanalytic.csv)
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(zoo)
library(tidyr)
library(lubridate)
library(tidyverse)
library(knitr)
library(kableExtra)
library(haven)
library(nimble)
library(stringi)
library(readxl)
library(openxlsx)
library(broom)


# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
# NOTE: Update the path below to match your local data directory.
# Example: combined_hit_analytic <- fread("data/hitanalytic.csv")

combined_hit_analytic <- fread("path/to/hitanalytic.csv")


# -----------------------------------------------------------------------------
# Data Cleaning & Variable Recoding
# -----------------------------------------------------------------------------

# Collapse Hispanic and missing race into a single "other/missing" category
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(
    patient_race_category = if_else(patient_race_category == "hispanic", "other", patient_race_category),
    patient_race_category = if_else(patient_race_category %in% c("other", "missing"), "other/missing", patient_race_category)
  )

# Replace missing/blank anticoagulation agent fields with "None"
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(
    daybeforehitlabacagent = if_else(is.na(daybeforehitlabacagent) | daybeforehitlabacagent == "", "None", daybeforehitlabacagent),
    dayafterhitlabacagent  = if_else(is.na(dayafterhitlabacagent)  | dayafterhitlabacagent  == "", "None", dayafterhitlabacagent)
  )

# Standardize gender labels
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(patient_gender = case_when(
    patient_gender %in% c("M", "MALE", "Male")     ~ "Male",
    patient_gender %in% c("F", "FEMALE", "Female") ~ "Female",
    TRUE ~ NA_character_
  ))

# Recode creatinine category, separating hemodialysis (HD) patients
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(creatinine_cat = case_when(
    HD == 1                                                              ~ "HD",
    creatinine_cat == "2.0+" & (is.na(HD) | HD == 0)                    ~ "Cr > 2, NOT on HD",
    creatinine_cat %in% c("<1.2", "1.2-2.0", "<2.0") & (is.na(HD) | HD == 0) ~ "Cr < 2",
    TRUE ~ NA_character_
  ))

table(combined_hit_analytic$creatinine_cat)

# Categorize HIT antibody lab result (numeric OD values and text results)
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(
    lab_result_clean  = str_trim(lab_result),
    lab_result_num    = suppressWarnings(as.numeric(lab_result_clean)),
    hit_ab_result_cat = case_when(
      # Numeric OD value categories
      !is.na(lab_result_num) & lab_result_num < 0.5                        ~ "<0.5 OR NEGATIVE",
      !is.na(lab_result_num) & lab_result_num >= 0.5 & lab_result_num < 1  ~ "0.5-1",
      !is.na(lab_result_num) & lab_result_num >= 1   & lab_result_num <= 1.8 ~ "1-1.8",
      !is.na(lab_result_num) & lab_result_num > 1.8                         ~ ">1.8",
      # String inequality notation (e.g., "<0.075")
      str_detect(lab_result_clean, "^<\\s*0\\.\\d+")                       ~ "<0.5 OR NEGATIVE",
      # Text-based results
      str_detect(tolower(lab_result_clean), "negative")                    ~ "<0.5 OR NEGATIVE",
      str_detect(tolower(lab_result_clean), "positive")                    ~ "positive WITHOUT numerical value",
      TRUE ~ NA_character_
    )
  )

table(combined_hit_analytic$hit_ab_result_cat)

combined_hit_analytic <- as_tibble(combined_hit_analytic)


# -----------------------------------------------------------------------------
# Baseline Characteristics Table
# -----------------------------------------------------------------------------

# Summary statistics for demographic and clinical variables
summary_table <- list(
  age = list(
    mean_sd = combined_hit_analytic %>%
      summarize(mean = mean(patient_age, na.rm = TRUE), sd = sd(patient_age, na.rm = TRUE)),
    range = range(combined_hit_analytic$patient_age, na.rm = TRUE)
  ),
  sex               = combined_hit_analytic %>% count(patient_gender)              %>% mutate(percent = n / sum(n) * 100),
  race              = combined_hit_analytic %>% count(patient_race_category)       %>% mutate(percent = n / sum(n) * 100),
  ethnicity         = combined_hit_analytic %>% count(patient_ethnicity_category)  %>% mutate(percent = n / sum(n) * 100),
  insurance         = combined_hit_analytic %>% count(primary_payer_type_category) %>% mutate(percent = n / sum(n) * 100),
  healthcare_system = combined_hit_analytic %>% count(source)                      %>% mutate(percent = n / sum(n) * 100),
  admission_service = combined_hit_analytic %>% count(admission_service_category)  %>% mutate(percent = n / sum(n) * 100)
)

remaining_stats <- list(
  platelet_count_day0          = combined_hit_analytic %>% count(platelet_count_cat)      %>% mutate(percent = n / sum(n) * 100),
  creatinine_day1              = combined_hit_analytic %>% count(creatinine_cat)           %>% mutate(percent = n / sum(n) * 100),
  anticoagulation_d1_intensity = combined_hit_analytic %>% count(daybeforehitlabaclvl)    %>% mutate(percent = n / sum(n) * 100),
  anticoagulation_d1_agent     = combined_hit_analytic %>% count(dayafterhitlabaclevel)   %>% mutate(percent = n / sum(n) * 100)
)

summary_table
remaining_stats


# -----------------------------------------------------------------------------
# Lab Result Turnaround Time (HIT Lab Result Lag)
# -----------------------------------------------------------------------------

# Calculate days from order to result
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(hit_result_lag_days = as.numeric(as.Date(hit_lab_result_date) - as.Date(hit_lab_date)))

# Overall turnaround time
summary_stats <- combined_hit_analytic %>%
  summarise(
    median   = median(hit_result_lag_days, na.rm = TRUE),
    IQR_low  = quantile(hit_result_lag_days, 0.25, na.rm = TRUE),
    IQR_high = quantile(hit_result_lag_days, 0.75, na.rm = TRUE),
    min      = min(hit_result_lag_days, na.rm = TRUE),
    max      = max(hit_result_lag_days, na.rm = TRUE)
  )
print(summary_stats)

# Turnaround time stratified by site
summary_stats_by_site <- combined_hit_analytic %>%
  group_by(source) %>%
  summarise(
    median   = median(hit_result_lag_days, na.rm = TRUE),
    IQR_low  = quantile(hit_result_lag_days, 0.25, na.rm = TRUE),
    IQR_high = quantile(hit_result_lag_days, 0.75, na.rm = TRUE),
    min      = min(hit_result_lag_days, na.rm = TRUE),
    max      = max(hit_result_lag_days, na.rm = TRUE),
    .groups  = "drop"
  )
print(summary_stats_by_site)


# -----------------------------------------------------------------------------
# Frequency Checks
# -----------------------------------------------------------------------------
table(combined_hit_analytic$HD)
table(combined_hit_analytic$daybeforehitlabaclvl)
table(combined_hit_analytic$daybeforehitlabacagent)
table(combined_hit_analytic$dayafterhitlabaclevel)
table(combined_hit_analytic$dayafterhitlabacagent)
table(combined_hit_analytic$accatdayprior)
table(combined_hit_analytic$creatinine_cat)
table(combined_hit_analytic$platelet_count_cat)
table(combined_hit_analytic$admission_cat_hit)
table(combined_hit_analytic$lab_result)

# Days from admission to HIT lab order (set negative values to NA)
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(days_after_admit_to_hitlab = ifelse(days_after_admit_to_hitlab < 0, NA, days_after_admit_to_hitlab))

print(median(combined_hit_analytic$days_after_admit_to_hitlab, na.rm = TRUE))

# NOTE: bundle_id may be an encounter-level identifier — verify de-identification
# before sharing any output that includes this variable.
ktableQ <- combined_hit_analytic %>%
  select(bundle_id, dayafterhitlabacagent, dayafterhitlabaclevel, low_dose_hep)


# =============================================================================
# Univariate Logistic Regression Models
# Outcome: outcome_dayafterhitlabacagent (anticoagulation agent day after HIT lab)
# =============================================================================

# Helper function to extract OR, 95% CI, and p-value from a glm model
extract_or <- function(model) {
  coef_mat <- summary(model)$coefficients
  data.frame(
    Variable  = rownames(coef_mat),
    OR        = round(exp(coef_mat[, "Estimate"]), 2),
    CI_Lower  = round(exp(coef_mat[, "Estimate"] - 1.96 * coef_mat[, "Std. Error"]), 2),
    CI_Upper  = round(exp(coef_mat[, "Estimate"] + 1.96 * coef_mat[, "Std. Error"]), 2),
    P_Value   = signif(coef_mat[, "Pr(>|z|)"], 3)
  )
}


# --- Age (continuous) --------------------------------------------------------
age_model <- glm(outcome_dayafterhitlabacagent ~ patient_age, data = combined_hit_analytic, family = binomial)
print(extract_or(age_model))


# --- Age (categorical: <65 vs 65+) -------------------------------------------
combined_hit_analytic$age_cat <- factor(
  ifelse(combined_hit_analytic$patient_age < 65, "<65", "65+")
)
combined_hit_analytic$age_cat <- relevel(combined_hit_analytic$age_cat, ref = "<65")

agecat_model <- glm(outcome_dayafterhitlabacagent ~ age_cat, data = combined_hit_analytic, family = binomial)
print(extract_or(agecat_model))


# --- Gender ------------------------------------------------------------------
combined_hit_analytic$patient_gender <- factor(combined_hit_analytic$patient_gender)
combined_hit_analytic$patient_gender <- relevel(combined_hit_analytic$patient_gender, ref = "Male")

gender_model <- glm(outcome_dayafterhitlabacagent ~ patient_gender, data = combined_hit_analytic, family = binomial)
print(extract_or(gender_model))


# --- Site --------------------------------------------------------------------
# Sites are anonymized; update reference level as appropriate.
combined_hit_analytic$site <- factor(combined_hit_analytic$source)
table(combined_hit_analytic$site)
combined_hit_analytic$site <- relevel(combined_hit_analytic$site, ref = "Site_1")

site_model <- glm(outcome_dayafterhitlabacagent ~ site, data = combined_hit_analytic, family = binomial)
print(extract_or(site_model))


# --- Admission Service -------------------------------------------------------
combined_hit_analytic$admission_cat_hit <- factor(combined_hit_analytic$admission_cat_hit)
combined_hit_analytic$admission_cat_hit <- relevel(combined_hit_analytic$admission_cat_hit, ref = "Non ICU Admission")

adm_model <- glm(outcome_dayafterhitlabacagent ~ admission_cat_hit, data = combined_hit_analytic, family = binomial)
print(extract_or(adm_model))


# --- Platelet Count ----------------------------------------------------------
combined_hit_analytic$platelet_count_cat <- factor(combined_hit_analytic$platelet_count_cat)
combined_hit_analytic$platelet_count_cat <- relevel(combined_hit_analytic$platelet_count_cat, ref = "150+")

plt_model <- glm(outcome_dayafterhitlabacagent ~ platelet_count_cat, data = combined_hit_analytic, family = binomial)
print(extract_or(plt_model))


# --- Creatinine --------------------------------------------------------------
combined_hit_analytic$creatinine_cat <- factor(combined_hit_analytic$creatinine_cat)
combined_hit_analytic$creatinine_cat <- relevel(combined_hit_analytic$creatinine_cat, ref = "<2.0")

cr_model <- glm(outcome_dayafterhitlabacagent ~ creatinine_cat, data = combined_hit_analytic, family = binomial)
print(extract_or(cr_model))


# --- Anticoagulation Intensity (Day Before HIT Lab) --------------------------
combined_hit_analytic$daybeforehitlabaclvl <- factor(combined_hit_analytic$daybeforehitlabaclvl)
combined_hit_analytic$daybeforehitlabaclvl <- relevel(combined_hit_analytic$daybeforehitlabaclvl, ref = "Prophylactic")

ac_model <- glm(outcome_dayafterhitlabacagent ~ daybeforehitlabaclvl, data = combined_hit_analytic, family = binomial)
print(extract_or(ac_model))


# --- Race --------------------------------------------------------------------
combined_hit_analytic$patient_race_category <- factor(combined_hit_analytic$patient_race_category)
combined_hit_analytic$patient_race_category <- relevel(combined_hit_analytic$patient_race_category, ref = "white")

race_model <- glm(outcome_dayafterhitlabacagent ~ patient_race_category, data = combined_hit_analytic, family = binomial)
print(extract_or(race_model))


# --- HIT Antibody Test Result ------------------------------------------------
combined_hit_analytic$hit_ab_result_cat <- factor(combined_hit_analytic$hit_ab_result_cat)
combined_hit_analytic$hit_ab_result_cat <- relevel(combined_hit_analytic$hit_ab_result_cat, ref = "0.5-1")

hit_model <- glm(outcome_dayafterhitlabacagent ~ hit_ab_result_cat, data = combined_hit_analytic, family = binomial)
print(extract_or(hit_model))
