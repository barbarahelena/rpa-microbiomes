## Compute 16S alpha diversity by ethnicity (throat and nose).
## 1. Load rarefied data, retain ethnicities N > 50, and calculate observed
##    richness, Shannon, and Simpson diversity.
## 2. Summarise by ethnicity; run Kruskal-Wallis, pairwise Wilcoxon,
##    and covariate-adjusted linear regression analyses.
## 3. Cache diversity values, summaries, and tests for the report and Figure 1.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(broom)
library(ggpubr)

## Setup
setwd(here::here())
dir.create("results/alpha_diversity", recursive = TRUE, showWarnings = FALSE)

## Keep only ethnicity groups with more than n=50 samples (matches Table 1)
keep_groups <- function(ps, min_n = 50) {
    counts <- table(sample_data(ps)$EthnicityTotal)
    names(counts)[counts > min_n]
}

## Covariates for linear regression
## Note: MigrationGen and ResidenceDuration_BA are excluded because they are
## structurally NA for all Dutch participants (migration-specific variables),
## which would drop the entire Dutch group from complete-case analysis.
covariates <- c("Age_FU", "Sex", "BMI_FU", "Smoking_FU", "Antibiotics_FU",
                "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU",
                "PM10_mean", "PM25_mean", "NO2_mean", "EC_mean", "Season",
                "SeqBatch")

## Helper: run linear regression for one metric
run_regression <- function(df, metric, covariates) {
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

## ---- Load data and calculate alpha diversity for each site ----
sites <- list(
    throat = readRDS("data/processed/ps_throat_rarefied.RDS"),
    nose   = readRDS("data/processed/ps_nose_rarefied.RDS")
)

for (site_name in names(sites)) {
    ps <- sites[[site_name]]

    ## Filter to ethnicity groups with N > 50 in this site
    ps <- subset_samples(ps, EthnicityTotal %in% keep_groups(ps))

    ## Compute alpha diversity metrics
    alpha_df <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson")) |>
        rownames_to_column("sample_id")

    ## Join with metadata
    meta <- sample_data(ps) |>
        as("data.frame") |>
        rownames_to_column("sample_id") |>
        select(sample_id, EthnicityTotal, all_of(covariates)) |>
        mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)))

    alpha_df <- alpha_df |>
        left_join(meta, by = "sample_id") |>
        filter(!is.na(EthnicityTotal))

    n_samples <- nrow(alpha_df)
    group_ns <- alpha_df |> count(EthnicityTotal, name = "n") |> arrange(desc(n))
    cat("Groups (N>50) for", site_name, ":", paste(group_ns$EthnicityTotal, collapse = ", "), "\n")

    ## Long format for plotting and summaries
    alpha_long <- alpha_df |>
        pivot_longer(cols = c(Observed, Shannon, Simpson),
                     names_to = "metric", values_to = "value")

    ## ---- Summary statistics by ethnicity ----
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

    ## ---- Kruskal-Wallis omnibus test across all groups ----
    kruskal_results <- alpha_long |>
        group_by(metric) |>
        summarise(
            statistic = kruskal.test(value ~ EthnicityTotal)$statistic,
            p.value   = kruskal.test(value ~ EthnicityTotal)$p.value,
            .groups   = "drop"
        )

    ## ---- Pairwise Wilcoxon rank-sum tests between every pair of groups ----
    pairwise_results <- alpha_long |>
        group_by(metric) |>
        group_modify(~ {
            pw <- pairwise.wilcox.test(.x$value, .x$EthnicityTotal, p.adjust.method = "BH")
            as.data.frame(as.table(pw$p.value)) |>
                filter(!is.na(Freq)) |>
                dplyr::rename(group1 = Var1, group2 = Var2, p.adj = Freq)
        }) |>
        ungroup()

    ## ---- Linear regression (adjusted for covariates) ----
    reg_results <- bind_rows(
        run_regression(alpha_df, "Observed", covariates),
        run_regression(alpha_df, "Shannon", covariates),
        run_regression(alpha_df, "Simpson", covariates)
    )
    ## Evaluate each facet's former ggpubr test during computation.
    kruskal_annotations <- alpha_long |>
        group_by(metric) |>
        group_modify(~ ggpubr::compare_means(value ~ EthnicityTotal,
                                            data = .x, method = "kruskal.test")) |>
        ungroup()
    ## Preserve Figure 1's original test calls and unrounded p-values.
    figure_tests <- list(
        kw_p = kruskal.test(Shannon ~ EthnicityTotal, data = alpha_df)$p.value,
        pw = pairwise.wilcox.test(alpha_df$Shannon, alpha_df$EthnicityTotal,
                                 p.adjust.method = "BH")
    )
    ## ---- Cache diversity values, summaries, and test results ----
    dir.create("results/alpha_diversity/cache", recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(alpha_df = alpha_df,
            alpha_long = alpha_long,
            n_samples = n_samples,
            summary_table = summary_table,
            reg_results = reg_results,
            group_ns = group_ns,
            kruskal_results = kruskal_results,
            pairwise_results = pairwise_results,
            kruskal_annotations = kruskal_annotations,
            figure_tests = figure_tests),
            paste0("results/alpha_diversity/cache/alpha_16s_", site_name, ".rds"))
}
