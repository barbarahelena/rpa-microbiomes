## Alpha diversity analysis: shotgun metagenomics (tongue and throat)
## Stratified by ethnicity with linear regression

## Libraries
library(here)
library(tidyverse)
library(vegan)
library(broom)

## Setup
setwd(here::here())
dir.create("results/alpha_diversity", recursive = TRUE, showWarnings = FALSE)

## Covariates for linear regression
## Note: MigrationGen and ResidenceDuration_BA are excluded because they are
## structurally NA for all Dutch participants (migration-specific variables),
## which would drop the entire Dutch group from complete-case analysis.
covariates <- c("Age_FU", "Sex", "BMI_FU", "Smoking_FU", "Antibiotics_FU",
                "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU",
                "PM10_mean", "PM25_mean", "NO2_mean", "EC_mean")

## Helper: compute alpha diversity from MetaPhlAn counts matrix (0-100 scale)
compute_alpha <- function(counts_mat) {
    ## Convert to proportions
    prop_mat <- counts_mat / 100
    ## Replace NA with 0
    prop_mat[is.na(prop_mat)] <- 0

    tibble(
        sample_id = rownames(counts_mat),
        Observed  = rowSums(prop_mat > 0),
        Shannon   = diversity(prop_mat, index = "shannon")
    )
}

## Helper: run linear regression for one metric
run_regression <- function(df, metric, covariates) {
    ## Restrict to complete cases and drop unused factor levels
    model_vars <- c(metric, "EthnicityTotal", covariates)
    df_cc <- df[complete.cases(df[, model_vars]), ] |>
        mutate(across(where(is.factor), droplevels))

    ## Drop covariates with < 2 unique values in complete-case subset
    usable <- covariates[sapply(covariates, function(v) {
        vals <- df_cc[[v]]
        if (is.factor(vals)) nlevels(vals) >= 2 else length(unique(vals)) >= 2
    })]
    dropped <- setdiff(covariates, usable)
    if (length(dropped) > 0)
        message("  Dropped (< 2 levels): ", paste(dropped, collapse = ", "))

    formula_str <- paste(metric, "~ EthnicityTotal +",
                         paste(usable, collapse = " + "))
    fit <- lm(as.formula(formula_str), data = df_cc)
    tidy(fit, conf.int = TRUE) |>
        mutate(metric = metric, .before = 1)
}

## ---- Analysis loop over sites ----
sites <- list(
    throat = readRDS("data/processed/shotgun_throat.RDS"),
    tongue = readRDS("data/processed/shotgun_tongue.RDS")
)

for (site_name in names(sites)) {
    site_data <- sites[[site_name]]

    ## Compute alpha diversity
    alpha_df <- compute_alpha(site_data$counts)

    ## Join with metadata
    meta <- site_data$sample_data |>
        rownames_to_column("sample_id") |>
        select(sample_id, EthnicityTotal, all_of(covariates))

    alpha_df <- alpha_df |>
        left_join(meta, by = "sample_id") |>
        filter(!is.na(EthnicityTotal))

    n_samples <- nrow(alpha_df)

    ## Long format for plotting
    alpha_long <- alpha_df |>
        pivot_longer(cols = c(Observed, Shannon),
                     names_to = "metric", values_to = "value")

    ## ---- Summary statistics ----
    summary_table <- alpha_long |>
        group_by(metric, EthnicityTotal) |>
        summarise(
            n      = n(),
            mean   = mean(value),
            sd     = sd(value),
            median = median(value),
            q25    = quantile(value, 0.25),
            q75    = quantile(value, 0.75),
            min    = min(value),
            max    = max(value),
            .groups = "drop"
        )

    ## ---- Linear regression ----
    reg_results <- bind_rows(
        run_regression(alpha_df, "Observed", covariates),
        run_regression(alpha_df, "Shannon", covariates)
    )
    dir.create("results/alpha_diversity/cache", recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(alpha_df = alpha_df,
            alpha_long = alpha_long,
            n_samples = n_samples,
            summary_table = summary_table,
            reg_results = reg_results),
            paste0("results/alpha_diversity/cache/alpha_shotgun_", site_name, ".rds"))
}
