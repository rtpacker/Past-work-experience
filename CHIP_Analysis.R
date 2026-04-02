# =============================================================================
# CHIP (Clonal Hematopoiesis of Indeterminate Potential) & Stroke Analysis
# Case-Cohort Design — REGARDS Cohort
# =============================================================================
# Description: Associations between CHIP mutations (2%, 5%, 10% VAF thresholds)
#              and incident stroke, using case-cohort Cox models (cchs design).
#
# Outcome:   stroke21 (incident stroke)
# Exposure:  CHIPCOMB2p, CHIPCOMB5p, CHIPCOMB10p
# Design:    Case-cohort; subcohort defined by random == "Y"
#

# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(haven)
library(survival)
library(Hmisc)
library(cchs)
library(car)


# -----------------------------------------------------------------------------
# Paths — update before running
# -----------------------------------------------------------------------------
data.path     <- "path/to/chip_data/"
function.path <- "path/to/common_functions/"

source(paste0(function.path, "hr.f_unitcontBycat.R"))
source(paste0(function.path, "sas.format.round.R"))


# =============================================================================
# Load Data
# =============================================================================
chipanalytic_sample_t <- fread(paste0(data.path, "CHIPANALYSISFINALDATA.csv"))


# =============================================================================
# VAF Tie Identification
# =============================================================================
# Flag IDs where more than one mutation shares rank_VAF_BL == 1

chipanalytic_filtered <- chipanalytic_sample_t %>%
  group_by(id) %>%
  mutate(vafblties = ifelse(sum(rank_VAF_BL == 1) > 1, 1, 0)) %>%
  ungroup()

chipanalytic_vaftie <- chipanalytic_filtered %>%
  filter(vafblties == 1, rank_VAF_BL == 1)

fwrite(chipanalytic_vaftie, paste0(data.path, "output/VAF_ties_rank1.csv"))


# =============================================================================
# Baseline Characteristics Table
# =============================================================================
table(chipanalytic_sample_t$random)

# One row per ID: stroke cases and random subcohort
chipstrokes <- chipanalytic_sample_t %>%
  filter(stroke21 == 1) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()

controls <- chipanalytic_sample_t %>%
  filter(random == "Y") %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()

table(chipstrokes$CHIPCOMB2p)

#removed data summary
#################



stroke_summary  <- summarize_data(chipstrokes)
control_summary <- summarize_data(controls)

results <- bind_rows(
  stroke_summary  %>% mutate(Group = "Stroke"),
  control_summary %>% mutate(Group = "Controls")
)
print(results)

table(controls$BMI_Cat)
table(chipstrokes$BMI_Cat)

fwrite(stroke_summary,  paste0(data.path, "output/stroke_summary_baseline.csv"))
fwrite(control_summary, paste0(data.path, "output/control_summary_baseline.csv"))


# =============================================================================
# Case-Cohort Setup (cchs)
# =============================================================================

# Factor coding with reference levels
chipanalytic_sample <- fread(paste0(data.path, "CHIPANALYSISFINALDATA.csv"))
chipanalytic_sample$Race        <- relevel(factor(chipanalytic_sample$Race,        levels = c("B", "W")), ref = "W")
chipanalytic_sample$CHIPCOMB2p  <- relevel(factor(chipanalytic_sample$CHIPCOMB2p,  levels = c("0", "1")), ref = "0")
chipanalytic_sample$CHIPCOMB5p  <- relevel(factor(chipanalytic_sample$CHIPCOMB5p,  levels = c("0", "1")), ref = "0")
chipanalytic_sample$CHIPCOMB10p <- relevel(factor(chipanalytic_sample$CHIPCOMB10p, levels = c("0", "1")), ref = "0")
chipanalytic_sample$age_cat65   <- as.factor(chipanalytic_sample$age_cat65)
chipanalytic_sample$Gender      <- as.factor(chipanalytic_sample$Gender)
chipanalytic_sample <- as_tibble(chipanalytic_sample)

# Sampling fractions: subcohort proportion within each stratum (wt)
frac_table <- chipanalytic_sample %>%
  group_by(wt) %>%
  summarize(
    n_total = n(),
    n_sub   = sum(random == "Y", na.rm = TRUE),
    frac    = n_sub / n_total,
    .groups = "drop"
  )
print(frac_table)

sampling_fracs        <- setNames(frac_table$frac, frac_table$wt)
chipanalytic_sample_t <- merge(chipanalytic_sample, frac_table %>% select(wt, frac),
                               by = "wt", all.x = TRUE)

table(chipanalytic_sample_t$CC_WEIGHT, chipanalytic_sample_t$random)

chipcomb2p_counts <- chipanalytic_sample %>%
  filter(random == "Y") %>%
  count(CHIPCOMB2p)


# =============================================================================
# CHIP 2% VAF Threshold — cchs Models
# =============================================================================

# Model 1: Demographics
model1_cchs_2p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race + CHIPCOMB2p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model1_cchs_2p)

# Model 2: + Cardiovascular risk factors
model2_cchs_2p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever + CHIPCOMB2p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model2_cchs_2p)

# Model 3: + Inflammatory and renal biomarkers
model3_cchs_2p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever +
    logcrp_std + EGFR_CKDEPI + CHIPCOMB2p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model3_cchs_2p)


# =============================================================================
# CHIP 5% VAF Threshold — cchs Models
# =============================================================================

# Model 1: Demographics
model1_cchs_5p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race + CHIPCOMB5p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model1_cchs_5p)

# Model 2: + Cardiovascular risk factors
model2_cchs_5p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever + CHIPCOMB5p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model2_cchs_5p)

# Model 3: + Inflammatory and renal biomarkers
model3_cchs_5p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever +
    logcrp_std + EGFR_CKDEPI + CHIPCOMB5p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model3_cchs_5p)


# =============================================================================
# CHIP 10% VAF Threshold — cchs Models
# =============================================================================

# Model 1: Demographics
model1_cchs_10p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race + CHIPCOMB10p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model1_cchs_10p)

# Model 2: + Cardiovascular risk factors
model2_cchs_10p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever + CHIPCOMB10p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model2_cchs_10p)

# Model 3: + Inflammatory and renal biomarkers
model3_cchs_10p <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever +
    logcrp_std + EGFR_CKDEPI + CHIPCOMB10p,
  data              = chipanalytic_sample_t,
  inSubcohort       = chipanalytic_sample_t$random == "Y",
  stratum           = chipanalytic_sample_t$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(model3_cchs_10p)


# =============================================================================
# Race Interaction Models — CHIP 2% VAF Threshold
# Standard Cox for Wald interaction test; cchs for inference
# =============================================================================

# Model 1 — interaction test
coxm1 <- coxph(
  Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race + CHIPCOMB2p + CHIPCOMB2p:Race,
  data = chipanalytic_sample, ties = "breslow"
)
Anova(coxm1, type = 3, test.statistic = "Wald")  # interaction p-value

coxm1_cchs <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + REGION + Race + Age:Race + CHIPCOMB2p + Race:CHIPCOMB2p,
  data              = chipanalytic_sample,
  inSubcohort       = chipanalytic_sample$random == "Y",
  stratum           = chipanalytic_sample$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(coxm1_cchs)

# Model 2 — interaction test
coxm2 <- coxph(
  Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever + CHIPCOMB2p + CHIPCOMB2p:Race,
  data = chipanalytic_sample, ties = "breslow"
)
Anova(coxm2, type = 3, test.statistic = "Wald")  # interaction p-value

coxm2_cchs <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + REGION + Race + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever + CHIPCOMB2p:Race,
  data              = chipanalytic_sample,
  inSubcohort       = chipanalytic_sample$random == "Y",
  stratum           = chipanalytic_sample$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(coxm2_cchs)

# Model 3 — interaction test
coxm3 <- coxph(
  Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + Race + REGION + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever +
    logcrp_std + EGFR_CKDEPI + CHIPCOMB2p + CHIPCOMB2p:Race,
  data = chipanalytic_sample, ties = "breslow"
)
Anova(coxm3, type = 3, test.statistic = "Wald")  # interaction p-value

coxm3_cchs <- cchs(
  formula = Surv(start_time, end_time, event = stroke21) ~
    Age + Gender + REGION + Race + Age:Race +
    SBP_std + Diab_SRMed_glu + Smoke + PCVD + Afib_SR_ECG +
    LVH_12 + Hyper_Meds_SR_ever +
    logcrp_std + EGFR_CKDEPI + CHIPCOMB2p:Race,
  data              = chipanalytic_sample,
  inSubcohort       = chipanalytic_sample$random == "Y",
  stratum           = chipanalytic_sample$wt,
  samplingFractions = sampling_fracs,
  precision         = 1
)
summary(coxm3_cchs)
