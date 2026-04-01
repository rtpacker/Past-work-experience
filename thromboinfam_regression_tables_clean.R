# =============================================================================
# Residential Segregation & Thrombo-Inflammatory Biomarkers
# Regression Tables
# =============================================================================
# Description: Survey-weighted linear and logistic regression models examining
#              associations between residential segregation indices and
#              thrombo-inflammatory biomarkers in the REGARDS cohort.
#
# Exposures:   Dissimilarity index, Interaction index, Isolation index
#              (each standardized by SD and split into tertiles)
#
# Outcomes:    log_ddimer, log_crp, log_ifn, log_tnf, log_il6 (log-transformed),
#              eselectin, fix (untransformed),
#              IL-1b (tertile comparisons: T1 vs T2, T1 vs T3)
#
# Models:      Unadjusted + 11 sequentially adjusted models (age/gender,
#              income, education, neighborhood SES, physical activity, diet,
#              stroke/TIA, CAD, smoking, alcohol, diabetes)
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
library(openxlsx)
library(survey)
library(broom)
library(purrr)


# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
# NOTE: Update the path below to match your local data directory.
# Example: analyticsamplethr <- readRDS("data/analyticsamplethr.rds")

analyticsamplethr <- readRDS("path/to/analyticsamplethr.rds")


# -----------------------------------------------------------------------------
# Exclusions: Remove missing values for key variables
# -----------------------------------------------------------------------------
exclude_vars <- c(
  "rs_dissimilarityb", "rs_interactionb", "rs_isolationb",
  "log_ddimer", "log_crp", "log_ifn", "log_tnf", "log_il6",
  "eselectin", "fix", "il1",
  "age", "gender", "race"
)

for (v in exclude_vars) {
  analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[[v]]) & !is.nan(analyticsamplethr[[v]]), ]
}


# -----------------------------------------------------------------------------
# Variable Creation: Standardized indices and tertiles
# -----------------------------------------------------------------------------

# Standardize residential segregation indices (divide by SD)
analyticsamplethr$dissimilarity_index_std <- analyticsamplethr$rs_dissimilarityb / sd(analyticsamplethr$rs_dissimilarityb)
analyticsamplethr$interaction_index_std   <- analyticsamplethr$rs_interactionb   / sd(analyticsamplethr$rs_interactionb)
analyticsamplethr$isolation_index_std     <- analyticsamplethr$rs_isolationb     / sd(analyticsamplethr$rs_isolationb)

# Tertiles for each index (reference = Tertile 1)
make_tertile_factor <- function(x) {
  f <- factor(ntile(x, 3), levels = c("1", "2", "3"))
  relevel(f, ref = "1")
}

analyticsamplethr$dissimilarity_index_tertiles <- make_tertile_factor(analyticsamplethr$rs_dissimilarityb)
analyticsamplethr$interaction_index_tertiles   <- make_tertile_factor(analyticsamplethr$rs_interactionb)
analyticsamplethr$isolation_index_tertiles     <- make_tertile_factor(analyticsamplethr$rs_isolationb)

# IL-1b: floor values below 0.01 to 0.009, then create tertile subsets
analyticsamplethr$il1b2 <- ifelse(analyticsamplethr$il1 < 0.01, 0.009, analyticsamplethr$il1)
analyticsamplethr$il1b2_tertiles <- cut(
  analyticsamplethr$il1b2,
  breaks = quantile(analyticsamplethr$il1b2, probs = c(0, 1/3, 2/3, 1)),
  labels = c("Tertile 1", "Tertile 2", "Tertile 3")
)

# Tertile 1 vs 2 subset
il1b1v2dat <- subset(analyticsamplethr, il1b2_tertiles %in% c("Tertile 1", "Tertile 2"))
il1b1v2dat$highest_tertile <- ifelse(il1b1v2dat$il1b2_tertiles == "Tertile 2", 1, 0)

# Tertile 1 vs 3 subset
il1b1v3dat <- subset(analyticsamplethr, il1b2_tertiles %in% c("Tertile 1", "Tertile 3"))
il1b1v3dat$highest_tertile <- ifelse(il1b1v3dat$il1b2_tertiles == "Tertile 3", 1, 0)


# =============================================================================
# Helper Functions
# =============================================================================

# List of model names used throughout
model_names <- c(
  "Unadjusted",
  "Age and gender adjusted",
  "age, gender, income",
  "age, gender, education",
  "age, gender, nSES",
  "age, gender, physical activity",
  "age, gender, diet",
  "age, gender, stroke/TIA",
  "age, gender, coronary artery disease",
  "smoking indvidual",
  "alcohol consumption indvidual",
  "diabetes  indvidual"
)

# Exposures
exposures <- c("dissimilarity_index_std", "interaction_index_std", "isolation_index_std")

# -----------------------------------------------------------------------------
# fit_svyglm_models(): Fit all 12 adjusted models for one exposure-biomarker pair
# Returns a named list of svyglm summaries.
# -----------------------------------------------------------------------------
fit_svyglm_models <- function(biomarker, exposure, data, family = gaussian()) {
  dgn <- svydesign(~1, prob = ~prob, strata = ~race * gender, data = data)

  formulas <- list(
    "Unadjusted"                           = paste(biomarker, "~", exposure),
    "Age and gender adjusted"              = paste(biomarker, "~", exposure, "+ age + gender"),
    "age, gender, income"                  = paste(biomarker, "~", exposure, "+ age + gender + income"),
    "age, gender, education"               = paste(biomarker, "~", exposure, "+ age + gender + ed_cat"),
    "age, gender, nSES"                    = paste(biomarker, "~", exposure, "+ age + gender + nses_tract"),
    "age, gender, physical activity"       = paste(biomarker, "~", exposure, "+ age + gender + exercise_cat"),
    "age, gender, diet"                    = paste(biomarker, "~", exposure, "+ age + gender + diet7"),
    "age, gender, stroke/TIA"              = paste(biomarker, "~", exposure, "+ age + gender + tia_sr"),
    "age, gender, coronary artery disease" = paste(biomarker, "~", exposure, "+ age + gender + cad_sr_ecg"),
    "smoking indvidual"                    = paste(biomarker, "~", exposure, "+ smoke + age + gender"),
    "alcohol consumption indvidual"        = paste(biomarker, "~", exposure, "+ alc_niaaa + age + gender"),
    "diabetes  indvidual"                  = paste(biomarker, "~", exposure, "+ diab_srmed_glu + age + gender")
  )

  lapply(formulas, function(f) {
    summary(svyglm(as.formula(f), design = dgn, family = family, data = data))
  })
}

# -----------------------------------------------------------------------------
# extract_row_linear(): Extract one result row for untransformed (linear) biomarkers
# Returns beta, SD, SE, p-value, 95% CI on the original scale.
# -----------------------------------------------------------------------------
extract_row_linear <- function(model, exposure, biomarker_sd) {
  coef_mat    <- model$coefficients
  beta        <- coef_mat[exposure, "Estimate"]
  std_err     <- coef_mat[exposure, "Std. Error"]
  p_raw       <- coef_mat[exposure, "Pr(>|t|)"]
  p_fmt       <- ifelse(round(p_raw, 3) == 0.000, "<0.001", sprintf("%.3f", p_raw))
  moe         <- qnorm(0.975) * std_err

  data.frame(
    Beta      = sprintf("%.3f", beta),
    SD_biomarker = sprintf("%.3f", biomarker_sd),
    Std_Error = sprintf("%.3f", std_err),
    P_Value   = p_fmt,
    CI_95     = sprintf("(%.3f, %.3f)", beta - moe, beta + moe),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# extract_row_log(): Extract one result row for log-transformed biomarkers.
# Returns beta, SE, p-value, and reverse-transformed % change + CI.
# -----------------------------------------------------------------------------
extract_row_log <- function(model, exposure, biomarker_sd) {
  coef_mat   <- model$coefficients
  beta       <- coef_mat[exposure, "Estimate"]
  std_err    <- coef_mat[exposure, "Std. Error"]
  p_raw      <- coef_mat[exposure, "Pr(>|t|)"]
  p_fmt      <- ifelse(round(p_raw, 3) == 0.000, "<0.001", sprintf("%.3f", p_raw))

  sign_b     <- ifelse(beta >= 0, 1, -1)
  rev_t      <- sign_b * (exp(abs(beta)) - 1) * 100

  ci_stat    <- 1.96 * std_err
  b_abs      <- abs(beta)
  lcit       <- (exp(b_abs - ci_stat) - 1) * 100
  ucit       <- (exp(b_abs + ci_stat) - 1) * 100
  lcl        <- sign_b * lcit
  ucl        <- sign_b * ucit
  if (lcl > ucl) { temp <- lcl; lcl <- ucl; ucl <- temp }

  data.frame(
    Beta      = sprintf("%.3f", beta),
    SD_biomarker = sprintf("%.3f", biomarker_sd),
    Std_Error = sprintf("%.3f", std_err),
    P_Value   = p_fmt,
    Reverse_T = sprintf("%.2f%%", rev_t),
    CI_pct    = sprintf("(%0.2f%%, %0.2f%%)", lcl, ucl),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# extract_row_or(): Extract one result row for logistic models (IL-1b tertiles).
# Returns beta, OR, SE, p-value, OR 95% CI.
# -----------------------------------------------------------------------------
extract_row_or <- function(model, exposure) {
  coef_mat  <- model$coefficients
  beta      <- coef_mat[exposure, "Estimate"]
  std_err   <- coef_mat[exposure, "Std. Error"]
  p_raw     <- coef_mat[exposure, "Pr(>|t|)"]
  p_fmt     <- ifelse(round(p_raw, 3) == 0.000, "<0.001", sprintf("%.3f", p_raw))
  moe       <- qnorm(0.975) * std_err

  data.frame(
    Beta      = sprintf("%.3f", beta),
    OR        = sprintf("%.3f", exp(beta)),
    Std_Error = sprintf("%.3f", std_err),
    P_Value   = p_fmt,
    OR_CI_95  = sprintf("(%.3f, %.3f)", exp(beta - moe), exp(beta + moe)),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# build_results_wb(): Build a 3-sheet workbook (Dissimilarity/Isolation/Interaction)
# for a single biomarker, given its model results and extract function.
# -----------------------------------------------------------------------------
build_results_wb <- function(model_results_list, biomarker_key, exposures,
                              extract_fn, ..., output_path) {
  wb <- createWorkbook()
  sheet_map <- list(
    Dissimilarity_Results = exposures[1],
    Isolation_Results     = exposures[3],
    Interaction_Results   = exposures[2]
  )

  for (sheet_name in names(sheet_map)) {
    exposure   <- sheet_map[[sheet_name]]
    result_key <- paste0(exposure, "_", biomarker_key)
    models     <- model_results_list[[result_key]]

    rows <- lapply(model_names, function(mn) {
      cbind(Model = mn, extract_fn(models[[mn]], exposure, ...))
    })
    df <- do.call(rbind, rows)

    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, df, colNames = TRUE)
  }

  saveWorkbook(wb, output_path, overwrite = TRUE)
  message("Saved: ", output_path)
}


# =============================================================================
# SECTION 1: Fit Models — All Exposures x Continuous Biomarkers
# =============================================================================
biomarkers_continuous <- c("log_ddimer", "log_crp", "log_ifn", "log_tnf", "log_il6", "eselectin", "fix")

model_results <- list()

for (exposure in exposures) {
  for (biomarker in biomarkers_continuous) {
    key <- paste0(exposure, "_", biomarker)
    model_results[[key]] <- fit_svyglm_models(biomarker, exposure, analyticsamplethr)
  }
}


# =============================================================================
# SECTION 2: Export Results — Log-Transformed Biomarkers
# (Reverse-transformed % change + CI)
# =============================================================================
log_biomarkers <- list(
  log_ddimer = "ddimer",
  log_crp    = "crp",
  log_ifn    = "ifn",
  log_tnf    = "tnf",
  log_il6    = "il6"
)

for (bm_key in names(log_biomarkers)) {
  bm_sd <- sd(analyticsamplethr[[sub("log_", "", bm_key)]])  # SD of original scale var

  build_results_wb(
    model_results_list = model_results,
    biomarker_key      = bm_key,
    exposures          = exposures,
    extract_fn         = extract_row_log,
    biomarker_sd       = bm_sd,
    output_path        = paste0("output/", log_biomarkers[[bm_key]], "_results_weighted.xlsx")
  )
}


# =============================================================================
# SECTION 3: Export Results — Untransformed Biomarkers (eselectin, FIX)
# (Beta + linear 95% CI)
# =============================================================================
linear_biomarkers <- c("eselectin", "fix")

for (bm_key in linear_biomarkers) {
  bm_sd <- sd(analyticsamplethr[[bm_key]])

  build_results_wb(
    model_results_list = model_results,
    biomarker_key      = bm_key,
    exposures          = exposures,
    extract_fn         = extract_row_linear,
    biomarker_sd       = bm_sd,
    output_path        = paste0("output/", bm_key, "_results_weighted.xlsx")
  )
}


# =============================================================================
# SECTION 4: IL-1b — Tertile Logistic Models (T1 vs T2, T1 vs T3)
# =============================================================================

# Helper to fit IL-1b logistic models for one tertile subset
fit_il1b_models <- function(data) {
  fit_svyglm_models("highest_tertile", exposures, data, family = binomial())
  # Note: fit_svyglm_models loops internally; wrap to iterate exposures
  results <- list()
  for (exp in exposures) {
    key <- paste0(exp, "_highest_tertile")
    results[[key]] <- fit_svyglm_models("highest_tertile", exp, data, family = binomial())
  }
  results
}

model_results_il1b1v2dat <- fit_il1b_models(il1b1v2dat)
model_results_il1b1v3dat <- fit_il1b_models(il1b1v3dat)

# Export IL-1b T1 vs T2
build_results_wb(
  model_results_list = model_results_il1b1v2dat,
  biomarker_key      = "highest_tertile",
  exposures          = exposures,
  extract_fn         = extract_row_or,
  output_path        = "output/il1b_t1vt2_results_weighted.xlsx"
)

# Export IL-1b T1 vs T3
build_results_wb(
  model_results_list = model_results_il1b1v3dat,
  biomarker_key      = "highest_tertile",
  exposures          = exposures,
  extract_fn         = extract_row_or,
  output_path        = "output/il1b_t1vt3_results_weighted.xlsx"
)
