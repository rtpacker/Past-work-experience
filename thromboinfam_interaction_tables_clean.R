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
# NOTE: Update data.path before running.
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(haven)
library(openxlsx)
library(survey)
library(broom)
library(purrr)


# -----------------------------------------------------------------------------
# Paths — update before running
# -----------------------------------------------------------------------------
data.path <- "path/to/data/"


# =============================================================================
# Load Data, Exclusions & Variable Derivation
# =============================================================================
analyticsamplethr <- readRDS(paste0(data.path, "analyticsamplethr.rds"))

# Remove missing values for key variables
exclude_vars <- c(
  "rs_dissimilarityb", "rs_interactionb", "rs_isolationb",
  "log_ddimer", "log_crp", "log_ifn", "log_tnf", "log_il6",
  "eselectin", "fix", "il1", "age", "gender", "race"
)
for (v in exclude_vars) {
  analyticsamplethr <- analyticsamplethr[
    !is.na(analyticsamplethr[[v]]) & !is.nan(analyticsamplethr[[v]]), ]
}

# Standardize residential segregation indices (divide by SD)
analyticsamplethr$dissimilarity_index_std <- analyticsamplethr$rs_dissimilarityb / sd(analyticsamplethr$rs_dissimilarityb)
analyticsamplethr$interaction_index_std   <- analyticsamplethr$rs_interactionb   / sd(analyticsamplethr$rs_interactionb)
analyticsamplethr$isolation_index_std     <- analyticsamplethr$rs_isolationb     / sd(analyticsamplethr$rs_isolationb)

# Tertiles for each index (reference = Tertile 1)
analyticsamplethr <- analyticsamplethr %>%
  mutate(
    dissimilarity_index_tertiles = ntile(rs_dissimilarityb, 3),
    interaction_index_tertiles   = ntile(rs_interactionb,   3),
    isolation_index_tertiles     = ntile(rs_isolationb,     3)
  )

make_tertile_factor <- function(x) relevel(factor(x, levels = c("1", "2", "3")), ref = "1")
analyticsamplethr$dissimilarity_index_tertiles <- make_tertile_factor(analyticsamplethr$dissimilarity_index_tertiles)
analyticsamplethr$interaction_index_tertiles   <- make_tertile_factor(analyticsamplethr$interaction_index_tertiles)
analyticsamplethr$isolation_index_tertiles     <- make_tertile_factor(analyticsamplethr$isolation_index_tertiles)

# Survey design: stratified by race x gender
dgn <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr)


# =============================================================================
# Interaction Models — Race, Gender, Race x Gender
# =============================================================================
# For each biomarker, three interaction models are fit using the interaction
# index as the primary exposure (age + gender adjusted):
#   1. Exposure x Race
#   2. Exposure x Gender
#   3. Exposure x Race x Gender
#
# Results are printed to console. The exposure variable below can be swapped
# to dissimilarity_index_std or isolation_index_std as needed.

biomarkers <- c("fix", "log_ddimer", "log_crp", "log_ifn", "log_tnf", "log_il6", "eselectin")

interaction_results <- list()

for (biomarker in biomarkers) {

  # Interaction with race
  int_race <- svyglm(
    as.formula(paste(biomarker, "~ interaction_index_std * race + age + gender")),
    design = dgn, data = analyticsamplethr
  )

  # Interaction with gender
  int_gender <- svyglm(
    as.formula(paste(biomarker, "~ interaction_index_std * gender + age + race")),
    design = dgn, data = analyticsamplethr
  )

  # Interaction with race and gender
  int_race_gender <- svyglm(
    as.formula(paste(biomarker, "~ interaction_index_std * race * gender + age")),
    design = dgn, data = analyticsamplethr
  )

  interaction_results[[biomarker]] <- list(
    race        = summary(int_race),
    gender      = summary(int_gender),
    race_gender = summary(int_race_gender)
  )

  message(sprintf("\n=== %s ===", biomarker))
  message("--- x Race ---");        print(summary(int_race))
  message("--- x Gender ---");      print(summary(int_gender))
  message("--- x Race x Gender ---"); print(summary(int_race_gender))
}
