# Table 1 HELIUS population

# Libraries
library(tidyverse)
library(tableone)

# Data
meta <- readRDS("data/processed/HELIUSmetadata_clean.RDS")
dim(meta)

# Shotgun subset
tonguemeta <- meta |> filter(!is.na(TongueSampleID)) |> droplevels()
throatmeta <- meta |> filter(!is.na(ThroatSampleID)) |> droplevels()
all(tonguemeta$ID %in% throatmeta$ID)

vars_table1 <- c(
  # Demografie / etniciteit
  "Sex",
  "Age_FU",
  "EthnicityTotal",
  "MigrationGen",
  "ResidenceDuration_BA",

  # Leefstijl
  "Smoking_FU",
  "AlcoholYN_FU",
  "BMI_FU",

  # Cardiometabool
  "SBP_FU",
  "DBP_FU",
  "HTSelfBP_FU",
  "DMSelfGluc_FU",
  "MetSyn_FU",

  # Medicatie (microbioom-relevant)
  "Antibiotics_FU",
  "Antihypertensiva_FU",
  "Lipidlowering_FU",

  # Psychosociaal (optioneel maar HELIUS-typisch)
  "DiscrMean_BA",

  # Mond- en neushygiëne (cruciaal voor orale/nasale swabs)
  "ToothBrushing_FU",
  "TongueBrushing_FU",
  "Mouthwash_FU",
  "OralHealth_FU",
  "Nasal_FU"
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
ids <- sample_names(ps)
str_length(ids[[1]])
str_length(meta$ID[1]) # these are different IDs which makes matching with 16S data impossible
