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
# NOTE: Update data.path and output.path before running.
#       Re-run the model-fitting loop and export block for each stratum
#       by swapping the active design object and subset (see comments below).
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
data.path   <- "path/to/data/"
output.path <- "path/to/output/stratification_tables/"


# =============================================================================
# Load Data, Exclusions & Variable Derivation
# =============================================================================
analyticsamplethr <- readRDS(paste0(data.path, "analyticsamplethr.rds"))

# Remove missing values for key variables
exclude_vars <- c(
  "rs_dissimilarityb", "rs_interactionb", "rs_isolationb",
  "log_ddimer", "log_crp", "log_ifn", "log_tnf", "log_il6",
  "eselectin", "fix", "il1", "age", "gender"
)
for (v in exclude_vars) {
  analyticsamplethr <- analyticsamplethr[
    !is.na(analyticsamplethr[[v]]) & !is.nan(analyticsamplethr[[v]]), ]
}

# Standardize residential segregation indices (divide by SD)
analyticsamplethr$dissimilarity_index_std <- analyticsamplethr$rs_dissimilarityb / sd(analyticsamplethr$rs_dissimilarityb)
analyticsamplethr$interaction_index_std   <- analyticsamplethr$rs_interactionb   / sd(analyticsamplethr$rs_interactionb)
analyticsamplethr$isolation_index_std     <- analyticsamplethr$rs_isolationb     / sd(analyticsamplethr$rs_isolationb)


# =============================================================================
# Stratified Subsets
# =============================================================================
analyticsamplethr_W  <- analyticsamplethr %>% filter(race   == "W")
analyticsamplethr_B  <- analyticsamplethr %>% filter(race   == "B")
analyticsamplethr_M  <- analyticsamplethr %>% filter(gender == "M")
analyticsamplethr_F  <- analyticsamplethr %>% filter(gender == "F")
analyticsamplethr_BM <- analyticsamplethr %>% filter(race == "B" & gender == "M")
analyticsamplethr_BF <- analyticsamplethr %>% filter(race == "B" & gender == "F")
analyticsamplethr_WM <- analyticsamplethr %>% filter(race == "W" & gender == "M")
analyticsamplethr_WF <- analyticsamplethr %>% filter(race == "W" & gender == "F")


# =============================================================================
# Survey Designs — One per Stratum
# =============================================================================
w_dgn          <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr_W)
b_dgn          <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr_B)
men_dgn        <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr_M)
women_dgn      <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr_F)
black_male_dgn <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr_BM)
black_fem_dgn  <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr_BF)
white_male_dgn <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr_WM)
white_fem_dgn  <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = analyticsamplethr_WF)


# =============================================================================
# Helper Functions
# =============================================================================

# -----------------------------------------------------------------------------
# extract_strat_row()
# Extracts beta, SD of biomarker, SE, p-value, reverse-transformed % change,
# and sign-corrected 95% CI from one svyglm summary for a given exposure.
# -----------------------------------------------------------------------------
extract_strat_row <- function(model, exposure, model_name, data, biomarker) {
  coef_mat   <- model$coefficients
  beta       <- coef_mat[exposure, "Estimate"]
  std_err    <- coef_mat[exposure, "Std. Error"]
  p_raw      <- coef_mat[exposure, "Pr(>|t|)"]
  p_fmt      <- ifelse(round(p_raw, 3) == 0.000, "<0.001", sprintf("%.3f", p_raw))

  sign_b     <- ifelse(beta >= 0, 1, -1)
  rev_t      <- sign_b * (exp(abs(beta)) - 1) * 100

  ci_lo_raw  <- abs(beta) - 1.96 * std_err
  ci_hi_raw  <- abs(beta) + 1.96 * std_err
  lcl        <- sign_b * (exp(ci_lo_raw) - 1) * 100
  ucl        <- sign_b * (exp(ci_hi_raw) - 1) * 100
  if (lcl > ucl) { temp <- lcl; lcl <- ucl; ucl <- temp }

  data.frame(
    Model     = model_name,
    SD        = round(sd(data[[biomarker]], na.rm = TRUE), 3),
    Beta      = sprintf("%.3f", beta),
    Std_Error = sprintf("%.3f", std_err),
    P_Value   = p_fmt,
    Reverse_T = sprintf("%.2f%%", rev_t),
    CI        = sprintf("(%0.2f%%, %0.2f%%)", lcl, ucl),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# export_strat_wb()
# Fits unadjusted and age-adjusted svyglm models for the specified biomarkers
# and exposures within a given stratum, then saves a 3-sheet Excel workbook.
#
# Args:
#   biomarkers   : character vector of outcome variable names
#   active_design: svydesign object for the stratum
#   active_data  : data frame subset for the stratum
#   filename     : output filename (saved to output.path)
#   age_only     : if TRUE, age-adjusted model uses only age (no gender);
#                  set TRUE for single-gender strata
# -----------------------------------------------------------------------------
export_strat_wb <- function(biomarkers, active_design, active_data,
                             filename, age_only = FALSE) {
  exposures <- c("dissimilarity_index_std", "interaction_index_std", "isolation_index_std")
  models_of_interest <- c("Unadjusted", "Age adjusted")

  exposure_map <- list(
    Dissimilarity_Results = "dissimilarity_index_std",
    Isolation_Results     = "isolation_index_std",
    Interaction_Results   = "interaction_index_std"
  )

  age_covar <- if (age_only) "+ age" else "+ age + gender"

  for (biomarker in biomarkers) {

    # Fit models for all exposures
    model_results <- list()
    for (exposure in exposures) {
      key <- paste0(exposure, "_", biomarker)
      model_results[[key]] <- list(
        Unadjusted    = summary(svyglm(
          as.formula(paste(biomarker, "~", exposure)),
          design = active_design, data = active_data
        )),
        `Age adjusted` = summary(svyglm(
          as.formula(paste(biomarker, "~", exposure, age_covar)),
          design = active_design, data = active_data
        ))
      )
    }

    # Build workbook
    wb <- createWorkbook()
    for (sheet_name in names(exposure_map)) {
      exposure   <- exposure_map[[sheet_name]]
      result_key <- paste0(exposure, "_", biomarker)
      models     <- model_results[[result_key]]

      rows <- lapply(models_of_interest, function(mn) {
        extract_strat_row(models[[mn]], exposure, mn, active_data, biomarker)
      })
      df <- do.call(rbind, rows)

      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, df, colNames = TRUE)
    }

    out_file <- paste0(output.path, sub("log_", "", biomarker), "_", filename, ".xlsx")
    saveWorkbook(wb, out_file, overwrite = TRUE)
    message("Saved: ", basename(out_file))
  }
}


# =============================================================================
# Export — One Call per Stratum
# =============================================================================
# Biomarkers of primary interest for stratified analyses
biomarkers_strat <- c("log_ddimer", "log_crp", "log_ifn", "log_tnf",
                      "log_il6", "eselectin", "fix")

# White participants
export_strat_wb(biomarkers_strat, w_dgn, analyticsamplethr_W,
                filename = "strat_white")

# Black participants
export_strat_wb(biomarkers_strat, b_dgn, analyticsamplethr_B,
                filename = "strat_black")

# Male participants (age_only = TRUE — no gender covariate in single-sex strata)
export_strat_wb(biomarkers_strat, men_dgn, analyticsamplethr_M,
                filename = "strat_male", age_only = TRUE)

# Female participants
export_strat_wb(biomarkers_strat, women_dgn, analyticsamplethr_F,
                filename = "strat_female", age_only = TRUE)

# Black Male
export_strat_wb(biomarkers_strat, black_male_dgn, analyticsamplethr_BM,
                filename = "strat_black_male", age_only = TRUE)

# Black Female
export_strat_wb(biomarkers_strat, black_fem_dgn, analyticsamplethr_BF,
                filename = "strat_black_female", age_only = TRUE)

# White Male
export_strat_wb(biomarkers_strat, white_male_dgn, analyticsamplethr_WM,
                filename = "strat_white_male", age_only = TRUE)

# White Female
export_strat_wb(biomarkers_strat, white_fem_dgn, analyticsamplethr_WF,
                filename = "strat_white_female", age_only = TRUE)
