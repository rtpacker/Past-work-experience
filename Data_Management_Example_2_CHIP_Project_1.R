# =============================================================================
# Data Management Example 2 — CHIP Project
# REGARDS Cohort: CHIP Mutation Data Pipeline
# =============================================================================
# Description: Full pipeline from raw sequencing data to the final analytic
#              dataset used for case-cohort Cox survival analysis.
#
# Pipeline overview:
#   Step 1 — Load reference datasets (REGARDS IDs, case-cohort weights)
#   Step 2 — K08 batch: ID linkage, VAF ranking
#   Step 3 — Firstset ASH batch: ID & timepoint linkage, VAF ranking
#   Step 4 — KL2 stroke batch: ID linkage, add CHIP-negative stroke participants
#   Step 5 — Combine all three sources into one CHIP analytic file
#   Step 6 — Merge covariates and stroke outcome
#   Step 7 — Derive analysis variables (CHIP thresholds, composite variables,
#             standardized biomarkers, age categories, factor coding)
#   Step 8 — Build case-cohort survival structure
#
# NOTE: Update data.path and ref.path to match your local directories.
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
library(readxl)


# -----------------------------------------------------------------------------
# Paths — update before running
# -----------------------------------------------------------------------------
data.path <- "path/to/chip_data/"   # CHIP project data directory
ref.path  <- "path/to/regards_ref/" # REGARDS reference/ID data directory


# =============================================================================
# STEP 1: Load Reference Datasets
# =============================================================================
# calcvars: main REGARDS covariate file (IDs, demographics, clinical variables)
# calcvars_ssurf: SSuRF follow-up visit covariates
# cc_labs_id: case-cohort design variables (random flag, CC_WEIGHT)

calcvars       <- read.sas7bdat(paste0(ref.path, "calcvars_id.sas7bdat"))
calcvars_ssurf <- read.sas7bdat(paste0(ref.path, "calcvars_ssurf.sas7bdat"))
cc_labs_id     <- read.sas7bdat(paste0(ref.path, "cc_labs_id.sas7bdat"))


# =============================================================================
# STEP 2: K08 Batch — ID Linkage & VAF Ranking
# =============================================================================
# The K08 batch contains the second set of sequenced stroke participants.
# All observations in this batch are baseline (BL) timepoint only.

# --- Per-person CHIP summary ---
k08_perp_raw <- read_xlsx(
  paste0(data.path, "Second set K08 remaining Stroke.xlsx"),
  sheet = "CH variants per person"
)
k08_perp <- k08_perp_raw %>%
  select(Sample, CH_incl_subclonal, N_CH_incl_subclonal,
         Top_CH_incl_subclonal_Gene, Top_CH_incl_subclonal_Gene_transcript,
         Top_CH_incl_subclonal_Gene_AA, Top_CH_incl_subclonal_Gene_VAF) %>%
  rename(Sample_ID = Sample)

# --- All-variants file ---
k08_av_raw <- read_xlsx(
  paste0(data.path, "Second set K08 remaining Stroke.xlsx"),
  sheet = "all_variants"
)
k08_av <- k08_av_raw %>%
  select(Sample, Gene.refGene, Chr, ExonicFunc.refGene, AAChange.refGene,
         cosmic70, transcriptOI, NonsynOI, VAF) %>%
  rename(Sample_ID = Sample)

# --- ID linkage ---
k08_link <- read_xlsx(paste0(data.path, "SecondSet_K08_remaining_stroke with ID linkage.xlsx")) %>%
  select(id, sample) %>%
  rename(Sample_ID = sample)

k08_perp <- merge(k08_perp, k08_link, by = "Sample_ID", all.x = TRUE) %>%
  as.data.table()

k08_av <- merge(k08_av, k08_link, by = "Sample_ID", all.x = TRUE) %>%
  as.data.table()

# Filter to IDs confirmed in REGARDS reference
k08_perp <- k08_perp[id %in% calcvars$id]
k08_av   <- k08_av[id %in% calcvars$id]

# Rank VAF within person (rank 1 = highest VAF; min method handles ties)
k08_av <- k08_av %>%
  group_by(id) %>%
  mutate(rank_VAF_BL   = rank(-VAF, na.last = "keep", ties.method = "min"),
         rank_VAF_SSuRF = NA_real_) %>%
  ungroup()

# Merge per-person summary onto all-variants
k08_allvars <- merge(k08_perp, k08_av, by = c("id", "Sample_ID"), all = TRUE)
k08_allvars$chip_assessed <- 1
k08_allvars$TIMEPOINT     <- "BL"
k08_allvars$DATA_source   <- "KO8"

fwrite(k08_perp, paste0(data.path, "output/k08_perp_final.csv"))
fwrite(k08_av,   paste0(data.path, "output/k08_av_final.csv"))


# =============================================================================
# STEP 3: Firstset ASH Batch — ID & Timepoint Linkage, VAF Ranking
# =============================================================================
# This batch contains both baseline (BL) and SSuRF follow-up observations.
# The "U" prefix on Sample IDs must be stripped to match the linkage file.

# --- ID + timepoint linkage ---
# ryan manifest provides packetid → timepoint mapping
ryan_manifest <- read_xlsx(
  paste0(data.path, "ryan first set Manifest IDs_pids_Timepoints 51023.xlsx")
) %>%
  rename(packetid = PacketID, id = Id_num) %>%
  select(id, packetid, TIMEPOINT)

id_link_raw <- read.sas7bdat(paste0(data.path, "ID_linkagefirst_crs_ash.sas7bdat"))
packetid_link <- read.sas7bdat(paste0(ref.path, "Idlink_packetid.sas7bdat")) %>%
  select(id, packetid)

ryan_manifest <- merge(ryan_manifest, packetid_link, by = "packetid", all.x = TRUE)

# Merge timepoint into ID linkage using id, Plate, Well as keys
ash_id_link <- merge(id_link_raw, ryan_manifest, by = c("id", "Plate", "Well")) %>%
  select(Sample_ID, id, TIMEPOINT)

# --- All-variants file ---
ash_av_raw <- read_xlsx(
  paste0(data.path, "Firstset_ASH CRS.xlsx"),
  sheet = "all_variants"
) %>%
  as.data.table()

ash_av <- ash_av_raw %>%
  select(Sample, Gene.refGene, Chr, ExonicFunc.refGene, cosmic70,
         transcriptOI, NonsynOI, VAF) %>%
  rename(Sample_ID = Sample)

# Strip leading "U" from Sample IDs to match linkage file
ash_av$Sample_ID <- sub("^U", "", ash_av$Sample_ID)

# Link to REGARDS IDs and timepoints
ash_av <- merge(ash_av, ash_id_link, by = "Sample_ID", all.x = TRUE)
table(ash_av$TIMEPOINT)  # BL and SSuRF counts

# Rank VAF separately within each timepoint
ash_av <- ash_av %>%
  group_by(id) %>%
  mutate(
    VAF_bl    = if_else(TIMEPOINT == "BL",    VAF, NA_real_),
    VAF_ssurf = if_else(TIMEPOINT == "SSuRF", VAF, NA_real_),
    rank_VAF_BL    = case_when(
      TIMEPOINT == "BL"    ~ rank(-VAF_bl,    na.last = "keep", ties.method = "min")
    ),
    rank_VAF_SSuRF = case_when(
      TIMEPOINT == "SSuRF" ~ rank(-VAF_ssurf, na.last = "keep", ties.method = "min")
    )
  ) %>%
  ungroup() %>%
  select(-VAF_bl, -VAF_ssurf)

# --- Per-person CHIP summary ---
ash_perp_raw <- read_xlsx(
  paste0(data.path, "Firstset_ASH CRS.xlsx"),
  sheet = "CHIP variants per person"
)
ash_perp <- ash_perp_raw %>%
  select(Sample, CHIP, N_CHIP, Top_CHIP_Gene, Top_CHIP_Gene_transcript,
         Top_CHIP_Gene_AA, Top_CHIP_Gene_VAF) %>%
  rename(Sample_ID = Sample)

ash_perp$Sample_ID <- sub("^U", "", ash_perp$Sample_ID)
ash_perp <- merge(ash_perp, ash_id_link, by = "Sample_ID", all.x = TRUE) %>%
  select(-Top_CHIP_Gene, -Top_CHIP_Gene_AA, -Top_CHIP_Gene_VAF, -Top_CHIP_Gene_transcript)

# Merge per-person onto all-variants
ash_allvars <- merge(ash_av, ash_perp, all = TRUE)
ash_allvars$chip_assessed <- 1
ash_allvars$DATA_source   <- "Firstset ASH"

fwrite(ash_allvars, paste0(data.path, "output/ash_allvars_final.csv"))


# =============================================================================
# STEP 4: KL2 Stroke Batch — ID Linkage & CHIP-Negative Participants
# =============================================================================
# KL2 contains CHIP-positive stroke participants. An additional stroke500
# file provides CHIP-negative stroke cases with no mutations detected.

# --- All-variants file ---
kl2_av_raw <- read_xlsx(paste0(data.path, "Stroke KL2.xlsx"))
kl2_link   <- read_xlsx(paste0(data.path, "KL2 Stroke Sample and ID linkage.xlsx"))

kl2_av <- kl2_av_raw %>%
  select(Sample, Gene.refGene, Chr, ExonicFunc.refGene, AAChange.refGene,
         cosmic70, transcriptOI, NonsynOI, VAF) %>%
  rename(`Sample ID` = Sample)

kl2_av <- merge(kl2_av, kl2_link, by = "Sample ID", all.x = TRUE) %>%
  as.data.table()

kl2_av$CHIP     <- 1
kl2_av$TIMEPOINT <- "BL"

# Rank VAF within person
kl2_av <- kl2_av %>%
  filter(TIMEPOINT == "BL") %>%
  group_by(id) %>%
  mutate(rank_VAF_BL = rank(-VAF, na.last = "keep", ties.method = "min")) %>%
  ungroup()

# --- Add CHIP-negative stroke participants from stroke500 linkage ---
# These are participants in the KL2 stroke sample with no mutations detected
stroke500 <- read.csv(paste0(data.path, "KL2 REGARDS 500 Stroke.csv"))

ids_chip_positive <- kl2_av$id
kl2_chip_neg <- stroke500 %>%
  filter(!(id %in% ids_chip_positive)) %>%
  select(id) %>%
  mutate(TIMEPOINT = "BL")

kl2_final <- rbind(kl2_av, kl2_chip_neg, fill = TRUE)
kl2_final <- kl2_final %>%
  rename(Sample_ID = `Sample ID`) %>%
  select(-`Packet ID`)

kl2_final$chip_assessed <- 1
kl2_final$DATA_source   <- "KL2"

fwrite(kl2_final, paste0(data.path, "output/kl2_final.csv"))


# =============================================================================
# STEP 5: Combine All Three Sources into One CHIP Analytic File
# =============================================================================
# Row-bind K08, Firstset ASH, and KL2 into a single long dataset.
# fill = TRUE handles columns that exist in some sources but not others.

chip_analytic <- rbindlist(
  list(ash_allvars, k08_allvars, kl2_final),
  fill = TRUE
)

fwrite(chip_analytic, paste0(data.path, "output/CHIPanalytic.csv"))

# Source counts
chip_analytic %>%
  group_by(DATA_source) %>%
  summarise(
    n_rows       = n(),
    n_unique_ids = n_distinct(id)
  )


# =============================================================================
# STEP 6: Merge Covariates and Stroke Outcome
# =============================================================================

# --- Case-cohort design variables: random flag ---
cc_vars <- cc_labs_id[cc_labs_id$id %in% chip_analytic$id, ] %>%
  dplyr::select(id, random)

chip_analytic <- merge(chip_analytic, cc_vars, by = "id", all.x = TRUE)

# Recode blank random values to "N"
chip_analytic$random[chip_analytic$random == ""] <- "N"
table(chip_analytic$random)

# --- Case-cohort weight (CC_WEIGHT) ---
fac9_fac11 <- read.sas7bdat(paste0(ref.path, "fixfxi_id.sas7bdat")) %>%
  select(id, CC_WEIGHT)
fac9_fac11 <- fac9_fac11[fac9_fac11$id %in% chip_analytic$id, ]

chip_analytic <- merge(chip_analytic, fac9_fac11, by = "id", all.x = TRUE)

# Verify subcohort count (expect ~800)
chip_analytic %>%
  group_by(random) %>%
  summarise(n_unique_ids = n_distinct(id))

# --- REGARDS covariate data ---
covar_vars <- c(
  "id", "Age", "Gender", "Race", "Insurance", "Income", "Income_4cat",
  "ED_Cat", "Ruca_4class", "Smoke", "Alc_Use", "PSS", "Diet7", "PA7",
  "Diabetes_SR", "Hyper_SR", "MI_SR_ECG", "BMI", "BMI_Cat", "Cholest",
  "EGFR_CKDEPI", "Crp", "Fram_CHD", "Reg_Nsaids", "Reg_Asa", "DVT_SR",
  "Death_indicator", "REGION", "Lipidemia_SR_Meds", "Hyper_SRmeds_BP",
  "Diab_SRMed_glu", "Afib_SR_ECG", "Cancer", "Wbc", "Pltc", "Hgb",
  "Rdwcv", "Mcv", "Albumin_serum", "Stroke_SR", "Hyper_Meds_SR_ever",
  "SBP", "PAD_amputation", "PAD_surgery", "LVH_12", "INTDATE",
  "Last_fudate", "InHomeDate"
)

calcvars_sub <- calcvars[calcvars$id %in% chip_analytic$id, covar_vars]
calcvars_sub$logcrp <- log(calcvars_sub$Crp)

chip_analytic <- merge(chip_analytic, calcvars_sub, by = "id", all.x = TRUE)

# --- Stroke outcome & follow-up time ---
# stroke22 provides the incident stroke indicator and event date
stroke22 <- read.csv(paste0(data.path, "stroke22_chipuvm.csv")) %>%
  filter(id %in% chip_analytic$id) %>%
  select(-Age, -Gender, -Race, -random)

chip_analytic <- merge(chip_analytic, stroke22, by = "id", all.x = TRUE)

# Convert SAS date integers to R Date format and compute follow-up time
chip_analytic <- chip_analytic %>%
  mutate(
    inhomedate  = as.Date(InHomeDate, origin = "1960-01-01"),
    intdate     = as.Date(INTDATE,    origin = "1960-01-01"),
    stroke21dt  = as.Date(stroke21dt),
    futime_stroke = as.numeric(difftime(stroke21dt, inhomedate, units = "days")) / 365.25
  )

fwrite(chip_analytic, paste0(data.path, "output/CHIPanalyticfinal.csv"))


# =============================================================================
# STEP 7: Derive Analysis Variables
# =============================================================================
chipanalytic_sample <- chip_analytic

# --- CHIP binary indicators at VAF thresholds ---
# VAF = 0 for participants with no detected mutation
chipanalytic_sample <- chipanalytic_sample %>%
  mutate(
    VAF        = ifelse(is.na(VAF), 0, VAF),
    CHIPCOMB2p  = ifelse(VAF >= 0.02, 1, 0),
    CHIPCOMB5p  = ifelse(VAF >= 0.05, 1, 0),
    CHIPCOMB10p = ifelse(VAF >= 0.10, 1, 0)
  )

# --- Combined CHIP indicator (ASH + K08 definitions) ---
chipanalytic_sample$CHIP[is.na(chipanalytic_sample$CHIP)] <- 0
chipanalytic_sample$CH_incl_subclonal[is.na(chipanalytic_sample$CH_incl_subclonal)] <- 0
chipanalytic_sample$CHIPCOMB <- ifelse(
  chipanalytic_sample$CHIP == 1 | chipanalytic_sample$CH_incl_subclonal == 1, 1, 0
)

# --- VAF rank: set NA to 1 (no competing mutation, so rank = 1) ---
chipanalytic_sample <- chipanalytic_sample %>%
  mutate(rank_VAF_BL = ifelse(is.na(rank_VAF_BL), 1, rank_VAF_BL))

# --- Peripheral vascular disease (PVD) composite ---
chipanalytic_sample$PVD <- NA
chipanalytic_sample$PVD[chipanalytic_sample$PAD_amputation == "Y" |
                          chipanalytic_sample$PAD_surgery   == "Y"] <- "Y"
chipanalytic_sample$PVD[chipanalytic_sample$PAD_amputation == "N" &
                          chipanalytic_sample$PAD_surgery   == "N"] <- "N"
chipanalytic_sample$PVD[chipanalytic_sample$PAD_amputation == ""  &
                          chipanalytic_sample$PAD_surgery   == ""] <- "N"
table(chipanalytic_sample$PVD)

# --- Prior cardiovascular disease: CHD or PVD ---
chipanalytic_sample$chd  <- chipanalytic_sample$MI_SR_ECG
chipanalytic_sample$PCVD <- ifelse(
  chipanalytic_sample$chd == "Y" | chipanalytic_sample$PVD == "Y", "Y", "N"
)
table(chipanalytic_sample$PCVD)

# --- Standardized biomarkers (divide by SD) ---
chipanalytic_sample$SBP_std    <- chipanalytic_sample$SBP    / sd(chipanalytic_sample$SBP,    na.rm = TRUE)
chipanalytic_sample$logcrp_std <- chipanalytic_sample$logcrp / sd(chipanalytic_sample$logcrp, na.rm = TRUE)

# --- Age categories ---
chipanalytic_sample$age_cat65 <- ifelse(chipanalytic_sample$Age < 65, "<65", "65+")
chipanalytic_sample$age_cat65 <- relevel(
  factor(chipanalytic_sample$age_cat65, levels = c("<65", "65+")), ref = "<65"
)

chipanalytic_sample <- chipanalytic_sample %>%
  mutate(age_cat_med = case_when(
    Age >= 45 & Age <= 50 ~ "45-50",
    Age >= 51 & Age <= 55 ~ "51-55",
    Age >= 56 & Age <= 60 ~ "56-60",
    Age >= 61 & Age <= 65 ~ "61-65",
    Age >= 66 & Age <= 70 ~ "66-70",
    Age >= 71 & Age <= 75 ~ "71-75",
    Age >= 76 & Age <= 80 ~ "76-80",
    Age  > 80             ~ ">80",
    TRUE ~ NA_character_
  ))

chipanalytic_sample$age_cat_med <- relevel(
  factor(chipanalytic_sample$age_cat_med,
         levels = c("45-50", "51-55", "56-60", "61-65", "66-70", "71-75", "76-80", ">80")),
  ref = "45-50"
)
table(chipanalytic_sample$age_cat_med)

# --- Factor coding with reference levels ---
chipanalytic_sample$Race       <- relevel(factor(chipanalytic_sample$Race,       levels = c("B", "W")), ref = "W")
chipanalytic_sample$CHIPCOMB    <- relevel(factor(chipanalytic_sample$CHIPCOMB,   levels = c("0", "1")), ref = "0")
chipanalytic_sample$CHIPCOMB2p  <- relevel(factor(chipanalytic_sample$CHIPCOMB2p,  levels = c("0", "1")), ref = "0")
chipanalytic_sample$CHIPCOMB5p  <- relevel(factor(chipanalytic_sample$CHIPCOMB5p,  levels = c("0", "1")), ref = "0")
chipanalytic_sample$CHIPCOMB10p <- relevel(factor(chipanalytic_sample$CHIPCOMB10p, levels = c("0", "1")), ref = "0")

chipanalytic_sample <- chipanalytic_sample %>%
  mutate(
    Gender = factor(Gender, levels = c("F", "M")),
    REGION = factor(REGION, levels = c("NONBELT", "BUCKLE", "BELT"))
  )


# =============================================================================
# STEP 8: Build Case-Cohort Survival Structure
# =============================================================================
# Case-cohort design requires two risk-interval rows for cases in the subcohort:
#   s1: one row per person covering their full follow-up
#       - cases:     enter at t-1 (just before event), weight = 1
#       - non-cases: enter at t = 0, weight = CC_WEIGHT
#   s2: additional row for cases also in the random subcohort
#       - covers t = 0 to t-1 with stroke21 = 0 (not yet a case), weight = CC_WEIGHT

anset <- chipanalytic_sample

# s1: one row per person
s1 <- anset %>%
  mutate(
    start_time = case_when(
      stroke21 == 1 ~ futime_stroke - 1,  # Cases enter just before event
      stroke21 == 0 ~ 0                   # Non-cases enter at time 0
    ),
    end_time = futime_stroke,
    wt = case_when(
      stroke21 == 1 ~ 1,          # Cases carry weight 1
      stroke21 == 0 ~ CC_WEIGHT   # Non-cases carry case-cohort weight
    )
  )

# s2: extra pre-event row for cases who are also in the random subcohort
s2_stroke <- anset %>%
  filter(stroke21 == 1 & random == "Y") %>%
  mutate(
    start_time = 0,
    end_time   = futime_stroke - 1,
    stroke21   = 0,         # Not yet a case during this interval
    wt         = CC_WEIGHT
  )

# Combine into final analytic dataset
data_frame_use <- bind_rows(s1, s2_stroke)

# Summary checks
data_frame_use %>% summarise(n_unique_ids = n_distinct(id))
data_frame_use %>% group_by(stroke21) %>% summarise(n = n(), n_ids = n_distinct(id))

fwrite(data_frame_use, paste0(data.path, "output/CHIPANALYSISFINALDATA.csv"))
