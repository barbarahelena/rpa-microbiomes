# Table 1 HELIUS population

# Libraries
library(here)
library(tidyverse)
library(tableone)
library(phyloseq)

# Set working directory to project root and create output directories if needed
setwd(here::here())
dir.create("results/tableone", recursive = TRUE, showWarnings = FALSE)

# Data
meta <- readRDS("data/processed/HELIUSmetadata_clean.RDS")
dim(meta)

# Save str() of every variable for reference (full variable overview/codebook)
str_lines <- unlist(lapply(names(meta), function(v) {
  c(paste0("== ", v, " =="), capture.output(str(meta[[v]])), "")
}))
writeLines(str_lines, "results/tableone/meta_str.txt")

# Keep only ethnicity groups with more than n=50 samples
drop_small_groups <- function(data, min_n = 50) {
  data |>
    add_count(EthnicityTotal, name = "n_group") |>
    filter(n_group > min_n) |>
    select(-n_group) |>
    droplevels()
}

# Shotgun subset
tonguemeta <- meta |> filter(!is.na(TongueSampleID)) |> droplevels() |> drop_small_groups()
throatmeta <- meta |> filter(!is.na(ThroatSampleID)) |> droplevels() |> drop_small_groups()
all(tonguemeta$ID %in% throatmeta$ID)

vars_table1 <- c(
  # Demographics
  "Sex",
  "Age_FU",
  "EthnicityTotal",
  "MigrationGen",
  "ResidenceDuration_BA",  

  # Cardiometabolic risk factors
  "Smoking_FU",
  "AlcoholYN_FU",
  "BMI_FU",
  "SBP_FU",
  "DBP_FU",
  "HTSelfBP_FU",
  "DMSelfGluc_FU",
  "MetSyn_FU",

  # Medication
  "Antibiotics_FU",
  "Antihypertensiva_FU",
  "Lipidlowering_FU",
  "Corticosteroids_FU",
  "SystemicSteroids_FU",
  "Antihistamines_FU",
  "DecongAllerg_FU",
  "Antidepressants_FU",
  "Psychotropics_FU", # any psychotropic med (SSRI/SNRI, TCA, antipsychotics, anxiolytics, hypnotics, mood stabilizers, stimulants, addiction meds) - relevant due to xerostomia/dry mouth effects on oral microbiome

  # Respiratory / allergy status (relevant to nose and throat microbiome)
  "AsthmaCOPD_FU",

  # Perceived Ethnic Discrimination score (only available at baseline)
  "DiscrMean_BA",

  # Mouth and nose variables
  "ToothBrushing_FU",
  "TongueBrushing_FU",
  "Mouthwash_FU",
  "OralHealth_FU",
  "Nasal_FU", # nasal medication
  "OwnTeeth_FU",
  "DentistVisit_FU",
  "DentistCavities_FU",
  "DentistCrown_FU",
  "DentistCheckup_FU",
  "DentistTartar_FU",
  "DentistInflammation_FU",
  "DentistOther_FU",

  # COVID-19 test results
  "Cov1_ResultText",
  "Cov2_ResultText"
)

table_one <- CreateTableOne(
  vars = vars_table1,
  strata = "EthnicityTotal",
  data = tonguemeta, # or throatmeta, doesn't matter
  test = TRUE
)
table_one_csv <- print(table_one, nonnormal = c(""), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table_one_csv, "results/tableone/table_shotgun.csv", row.names = TRUE)

## 16S subset
ps <- readRDS("data/processed/ps_throat_rarefied.RDS")
psmeta <- as(sample_data(ps), "data.frame") |> drop_small_groups()

table_one_16s <- CreateTableOne(
  vars = vars_table1,
  strata = "EthnicityTotal",
  data = psmeta,
  test = TRUE
)
table_one_16s_csv <- print(table_one_16s, nonnormal = c(""), quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing = TRUE)
write.csv(table_one_16s_csv, "results/tableone/table_16S.csv", row.names = TRUE)
