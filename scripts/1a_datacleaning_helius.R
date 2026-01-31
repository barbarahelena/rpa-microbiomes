## Data cleaning of clinical data

## Libraries
library(tidyverse)

# Change type of variable - automatic detection of capitalization
yesno_auto <- function(x) {
  vals <- as.character(x)
  if (any(c("Nee", "Ja") %in% vals, na.rm = TRUE)) {
    fct_recode(as.factor(x), "No" = "Nee", "Yes" = "Ja")
  } else if (any(c("nee", "ja") %in% vals, na.rm = TRUE)) {
    fct_recode(as.factor(x), "No" = "nee", "Yes" = "ja")
  } else {
    as.factor(x)
  }
}

zero_one <- function(x) fct_recode(x, "No"="0", "Yes"="1")

## HELIUS data
meta <- haven::read_sav("data/250606_HELIUS data Barbara Verhaar.sav")
names(meta)

df_new <- meta |> 
    dplyr::select(
                  # IDs and sample identifiers
                  ID, TongueSampleID, ThroatSampleID,
                  
                  # H1 Baseline demographics and questionnaires
                  Questionnaire_BA = H1_VragenlijstJN,
                  PhysicalExam_BA = H1_LichamelijkOnderzoekJN,
                  DatePhysExam_BA = H1_DatumLO,
                  Sex = H1_geslacht,
                  Age_BA = H1_lft,
                  MigrationGen = H1_MigrGeneratie,
                  Ethnicity = H1_etniciteit,
                  EthnicityTotal = H1_EtnTotaal,
                  
                  # H1 Discrimination
                  DiscrSum_BA = H1_Discr_sumscore,
                  DiscrMean_BA = H1_Discr_meanscore,
                  AnyDiscrim_BA = H1_AnyDiscrim,
                  
                  # H1 Smoking and substances
                  Fagerstrom_BA = H1_Fagerstrom,
                  FagerstromCat_BA = H1_FagerstromCat,
                  HSI_BA = H1_HSI,
                  AuditAlcohol_BA = H1_Auditalcohol,
                  AlcoholProblem_BA = H1_AlcProbleem,
                  AlcoholFunction_BA = H1_AlcFunctioneren,
                  AlcoholPeriods_BA = H1_AlcPeriodes,
                  AlcoholPeriodsAge_BA = H1_AlcPeriodesLft,
                  CuditCannabis_BA = H1_CuditCannabis,
                  
                  # H1 Wellbeing and depression
                  Depressed_BA = H1_WlbvSomber,
                  NoPleasure_BA = H1_WlbvGeenPlezier,
                  Impaired_BA = H1_WlbvBelemmerd,
                  DepAgePeriod_BA = H1_WlbvLftdPeriode,
                  DepFreqPeriod_BA = H1_WlbvFreqPeriode,
                  
                  # H1 Social support
                  SSQT_BA = H1_SSQT,
                  SSQSa_BA = H1_SSQSa,
                  SSQSb_BA = H1_SSQSb,
                  
                  # H1 Migration and acculturation
                  ResidenceDuration_BA = H1_ResidenceDuration,
                  AgeMigration_BA = H1_AgeMigration,
                  DifficultyDutch_BA = H1_MoeiteNLtaal,
                  CultFeelBerrys_BA = H1_CultVoelenBerrysAcc,
                  CultOrientBerrys_BA = H1_CultOrientBerrysAcc,
                  CultNetworkBerrys_BA = H1_CultNetworkBerrysAcc,
                  CultDistMeanScore0_BA = H1_CultDistMeanScore0_corrected,
                  CultDistMeanScore6_BA = H1_CultDistMeanScore6_corrected,
                  
                  # H2 Follow-up basics
                  PhysicalExam_FU = H2_LichamelijkOnderzoekJN,
                  DatePhysExam_FU = H2_DatumLO,
                  Age_FU = H2_lft,
                  FollowUpTime = H2_fu_time,
                  ResponseTotal_FU = H2_responstotaal,
                  MaritalStatus_FU = H2_Burgsta,
                  WorkParticipation_FU = H2_Arbeidsparticipatie,
                  
                  # H2 Smoking and alcohol
                  Smoking_FU = H2_Roken_H1combined,
                  AlcoholYN_FU = H2_AlcoholJN,
                  AlcoholConsumption_FU = H2_AlcoholConsumption,
                  
                  # H2 Physical measurements
                  BMI_FU = H2_LO_BMI,
                  WHR_FU = H2_LO_WHR,
                  SBP_FU = H2_LO_GemBPSysZit,
                  DBP_FU = H2_LO_GemBPDiaZit,
                  HR_FU = H2_LO_GemBPHRZit,
                  
                  # H2 Hypertension
                  HTSelf_FU = H2_HT_Self_H1combined,
                  HTBP_FU = H2_HT_BP,
                  HTSelfBP_FU = H2_HT_SelfBP_H1combined,
                  HTSelfBPMed_FU = H2_HT_SelfBPMed_H1combined,
                  HTBPMed_FU = H2_HT_BPMed,
                  Antihypertensiva_FU = H2_Antihypertensiva,
                  AntihypertensivaC02_FU = H2_AntihypertensivaC02,
                  AntihypertensivaC03_FU = H2_AntihypertensivaC03,
                  AntihypertensivaC07_FU = H2_AntihypertensivaC07,
                  AntihypertensivaC08_FU = H2_AntihypertensivaC08,
                  AntihypertensivaC09_FU = H2_AntihypertensivaC09,
                  
                  # H2 Diabetes
                  DMSelf_FU = H2_DM_Self_H1combined,
                  DMGluc_FU = H2_DM_Gluc,
                  DMSelfGluc_FU = H2_DM_SelfGluc_H1combined,
                  DMGlucMed_FU = H2_DM_GlucMed,
                  DMSelfGlucMed_FU = H2_DM_SelfGlucMed_H1combined,
                  DiabetesMeds_FU = H2_Diabetesmiddelen,
                  
                  # H2 Cardiovascular events
                  StrokeLoss_FU = H2_UitvalBer_H1combined,
                  Infarction_FU = H2_Infarct_H1combined,
                  CardiacIntervention_FU = H2_OperDotByp_H1combined,
                  Lipidlowering_FU = H2_Antilipaemica,
                  
                  # H2 Metabolic syndrome
                  MetSynCentralObesity_FU = H2_MetSyn_CentralObesity,
                  MetSynHighTG_FU = H2_MetSyn_HighTriglyceride,
                  MetSynLowHDL_FU = H2_MetSyn_LowHDL,
                  MetSynHighBP_FU = H2_MetSyn_HighBP,
                  MetSynHighGluc_FU = H2_MetSyn_HighGluc,
                  MetSyn_FU = H2_MetSyn_MetabolicSyndrome,
                  
                  # H2 Kidney function
                  eCC_FU = H2_eCC,
                  MDRD_eGFR_FU = H2_MDRD_eGFR,
                  CKDEPI_eGFR_FU = H2_CKDEPI_eGFR,
                  CKDEPI_stage_FU = H2_CKDEPI_stage,
                  Microalbuminuria_FU = H2_Microalbuminurie,
                  ACR_KDIGO_FU = H2_ACR_KDIGO,
                  CKDrisk_KDIGO_FU = H2_CKDrisk_KDIGO,
                  CKDEPIcr21_eGFR_FU = H2_CKDEPIcr21_eGFR_corrected,
                  CKDEPIcr21_stage_FU = H2_CKDEPIcr21_stage_corrected,
                  CKDrisk_KDIGO_cr21_FU = H2_CKDrisk_KDIGO_cr21_corrected,
                  
                  # H2 Other medications
                  Antibiotics_FU = H2_Antibiotica,
                  Antithrombotics_FU = H2_Antithrombotica,
                  Corticosteroids_FU = H2_Corticosteroiden,
                  DecongAllerg_FU = H2_DecongAllerg,
                  Nasal_FU = H2_Nasal,
                  AsthmaCOPD_FU = H2_AstmCOPD,
                  Antihistamines_FU = H2_Antihistaminica,
                  SystemicSteroids_FU = H2_SystSteroiden,
                  
                  # H2 CVD risk scores
                  SCORE1_CVDmort_FU = H2_SCORE1_RM_CVDmort,
                  SCORE1_CVDmort_NL_FU = H2_SCORE1_RM_CVDmort_NL,
                  SCORE1_CVDtot_NL_FU = H2_SCORE1_RM_CVDtot_NL,
                  SCORE1_CVDmort_TC_FU = H2_SCORE1_RM_CVDmort_TC,
                  SCORE1_CVDmort_TC_NL_FU = H2_SCORE1_RM_CVDmort_TC_NL,
                  SCORE1_CVDtot_TC_NL_FU = H2_SCORE1_RM_CVDtot_TC_NL,
                  SCORE2_CVDmort_UnCal_FU = H2_SCORE2_RM_CVDmort_UnCal,
                  SCORE2_CVDmort_Cal_FU = H2_SCORE2_RM_CVDmort_Cal,
                  FramScore_FU = H2_FramScore,
                  Fram_CVD_FU = H2_Fram_CVD,
                  
                  # H2 Psychotropic medications
                  Antipsychotics_FU = H2_Antipsychotica,
                  Anxiolytics_FU = H2_Anxiolytica,
                  Hypnotics_FU = H2_Hypnotica,
                  SSRI_SNRI_FU = H2_SSRI_SNRI,
                  MoodStabilizers_FU = H2_Stemmingsstabilisatoren,
                  Stimulants_FU = H2_Stimulantia,
                  TCA_FU = H2_TCA,
                  AddictionMeds_FU = H2_Verslavingsmedicatie,
                  Psychotropics_FU = H2_Psychotroop,
                  Antidepressants_FU = H2_Antidepressiva,
                  
                  # H2 Depression scores
                  PHQ9_sum_FU = H2_PHQ9_sumscore,
                  PHQ9_deprsymp_FU = H2_PHQ9_deprsymp,
                  MDD_FU = H2_MDD,
                  
                  # H2 Lab results
                  Glucose_FU = H2_Lab_UitslagGLUC,
                  HbA1c_FU = H2_Lab_UitslagIH1C,
                  Triglycerides_FU = H2_Lab_UitslagTRIG,
                  Cholesterol_FU = H2_Lab_UitslagCHOL,
                  HDL_FU = H2_Lab_UitslagHDLS,
                  LDL_FU = H2_Lab_UitslagRLDL,
                  Creatinine_HP_FU = H2_Lab_UitslagKREA_HP,
                  Microalbumin_FU = H2_Lab_UitslagMIAL,
                  Creatinine_UP_FU = H2_Lab_UitslagKREA_UP,
                  MicroalbCreat_FU = H2_Lab_UitslagMIKR,
                  
                  # H2 Swabs and oral health
                  SwabsDate_FU = H2swabs_datum,
                  ToothBrushing_FU = H2_Tandenpoetsen,
                  TongueBrushing_FU = H2_Tongpoetsen,
                  Mouthwash_FU = H2_Mondspoeling,
                  OralHealth_FU = H2_MondGezondheid,
                  OwnTeeth_FU = H2_Eigengebit,
                  DentistVisit_FU = H2_TandartsBezoek,
                  DentistCavities_FU = H2_TandartsGaatjes,
                  DentistCrown_FU = H2_TandartsKroon,
                  DentistCheckup_FU = H2_TandartsControle,
                  DentistTartar_FU = H2_TandartsTandsteen,
                  DentistInflammation_FU = H2_TandartsOntsteking,
                  DentistOther_FU = H2_TandartsOverig,
                  
                  # COVID-19 tests
                  Cov1_Date = Cov1_Date,
                  Cov1_Age = Cov1_Age,
                  Cov1_ResultText = Cov1_UitslagTekst,
                  Cov1_ResultNum = Cov1_UitslagNum,
                  Cov2_Date = Cov2_Date,
                  Cov2_Age = Cov2_Age,
                  Cov2_ResultText = Cov2_UitslagTekst,
                  Cov2_ResultNum = Cov2_UitslagNum
    )

df_new2 <- df_new |> 
    mutate(
           # Standardize missing values
           across(where(is.character), \(x) {
             ifelse(x %in% c("Missing", "Missing: n.v.t.", "niet ingevuld", "nvt", 
                            "No lab result", "Missing: not applicable", "missing",
                            "Missing: not measured", 
                            "See comments lab results", "Low (<1 mmol/L)",
                            "Low (<0,08 mmol/L)", "Low (<0,10 mmol/L)",
                            "Smoking status unknown", "Number or duration unknown",
                            "Rookstatus onbekend", "Rookduur en/of aantal onbekend"),
                    NA_character_, x)
           }),
           
           # Create sample IDs
           fecessample_FU = str_c("HELIFU_", ID),
           ID = str_c("S", as.character(ID)),
           
           # Convert numeric variables
           across(c("Age_BA", "Age_FU", "FollowUpTime", 
                    "DiscrSum_BA", "DiscrMean_BA", "ResidenceDuration_BA", "AgeMigration_BA",
                    "CultDistMeanScore0_BA", "CultDistMeanScore6_BA",
                    "Fagerstrom_BA", "FagerstromCat_BA", "HSI_BA", "AuditAlcohol_BA", "CuditCannabis_BA",
                    "SSQT_BA", "SSQSa_BA", "SSQSb_BA",
                    "BMI_FU", "WHR_FU", "SBP_FU", "DBP_FU", "HR_FU",
                    "MDRD_eGFR_FU", "CKDEPI_eGFR_FU", "CKDEPIcr21_eGFR_FU",
                    "SCORE1_CVDmort_FU", "SCORE1_CVDmort_NL_FU", "SCORE1_CVDtot_NL_FU",
                    "SCORE1_CVDmort_TC_FU", "SCORE1_CVDmort_TC_NL_FU", "SCORE1_CVDtot_TC_NL_FU",
                    "SCORE2_CVDmort_UnCal_FU", "SCORE2_CVDmort_Cal_FU",
                    "FramScore_FU", "Fram_CVD_FU",
                    "PHQ9_sum_FU",
                    "Glucose_FU", "HbA1c_FU", "Triglycerides_FU", "Cholesterol_FU",
                    "HDL_FU", "LDL_FU", "Creatinine_HP_FU", "Microalbumin_FU",
                    "Creatinine_UP_FU", "MicroalbCreat_FU",
                    "Cov1_Age", "Cov2_Age"), as.numeric),
           
           # Convert labelled variables to factors
           across(where(haven::is.labelled), \(x) haven::as_factor(x, levels = "labels")),
           
           # Automatically recode yes/no variables (detects Nee/Ja or nee/ja)
           across(c("Questionnaire_BA", "PhysicalExam_BA", "PhysicalExam_FU",
                    "AnyDiscrim_BA", "Depressed_BA", "NoPleasure_BA", "Impaired_BA",
                    "HTSelf_FU", "HTBP_FU", "HTSelfBP_FU", "HTSelfBPMed_FU", "HTBPMed_FU",
                    "Antihypertensiva_FU", "AntihypertensivaC02_FU", "AntihypertensivaC03_FU",
                    "AntihypertensivaC07_FU", "AntihypertensivaC08_FU", "AntihypertensivaC09_FU",
                    "DMSelf_FU", "DMGluc_FU", "DMSelfGluc_FU", "DMGlucMed_FU", "DMSelfGlucMed_FU",
                    "DiabetesMeds_FU", "StrokeLoss_FU", "Infarction_FU", "CardiacIntervention_FU",
                    "Lipidlowering_FU", "MetSynCentralObesity_FU", "MetSynHighTG_FU",
                    "MetSynLowHDL_FU", "MetSynHighBP_FU", "MetSynHighGluc_FU", "MetSyn_FU",
                    "Microalbuminuria_FU", "Antibiotics_FU", "Antithrombotics_FU",
                    "Corticosteroids_FU", "Antipsychotics_FU", "Anxiolytics_FU",
                    "Hypnotics_FU", "SSRI_SNRI_FU", "MoodStabilizers_FU", "Stimulants_FU",
                    "TCA_FU", "AddictionMeds_FU", "Psychotropics_FU", "Antidepressants_FU",
                    "AlcoholYN_FU"), yesno_auto),
           
           # Recode Ethnicity
           EthnicityTotal = forcats::fct_recode(
             EthnicityTotal, 
             "Other" = "Other/unknown",
             "Other" = "Other/unknown Surinamese"
           ),
           
           # Recode Smoking
           Smoking_FU = fct_recode(
             Smoking_FU, 
             "Former smoking" = "No, but I used to smoke",
             "Never" = "No, I have nevver smoked"
           ),
           
           # Recode Alcohol Consumption
           AlcoholConsumption_FU = fct_recode(
             AlcoholConsumption_FU,
             "low (men 0-4 gl/w, women 0-2 gl/w)" = "Low (men 0-4 gl/w, women 0-2 gl/w)",
             "moderate (men 5-14 gl/w, women 3-7 gl/w)" = "Moderate (men 5-14 gl/w, women 3-7 gl/w)",
             "high (men >14 gl/w, women >7 gl/wk)" = "High (men >14 gl/w, women >7 gl/wk)"
           )
    ) |> 
    droplevels()

dim(df_new2)

saveRDS(df_new2, "data/HELIUSmetadata_clean.RDS")
  