## Alpha diversity analysis: 16S microbiome (throat and nose)
## Stratified by ethnicity (Dutch vs South-Asian Surinamese)
## with Wilcoxon tests and covariate-adjusted linear regression

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(broom)
library(ggpubr)

## Functions
theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

## Setup
setwd(here::here())
dir.create("results/alpha_diversity", recursive = TRUE, showWarnings = FALSE)

## Define ethnicity colours and covariates
eth_colours <- c("Dutch" = "#1F78B4", "South-Asian Surinamese" = "#E31A1C")

## Covariates for linear regression
## Note: MigrationGen and ResidenceDuration_BA are excluded because they are
## structurally NA for all Dutch participants (migration-specific variables),
## which would drop the entire Dutch group from complete-case analysis.
covariates <- c("Age_FU", "Sex", "BMI_FU", "Smoking_FU", "Antibiotics_FU",
                "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU")

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

## ---- Analysis loop over sites ----
sites <- list(
    throat = readRDS("data/processed/ps_throat_rarefied.RDS"),
    nose   = readRDS("data/processed/ps_nose_rarefied.RDS")
)

rarefaction_depths <- c(throat = 7500, nose = 2000)

for (site_name in names(sites)) {
    ps <- sites[[site_name]]
    rare_depth <- rarefaction_depths[[site_name]]

    ## Filter to Dutch and South-Asian Surinamese only
    ps <- subset_samples(ps, EthnicityTotal %in% c("Dutch", "South-Asian Surinamese"))

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
    n_dutch <- sum(alpha_df$EthnicityTotal == "Dutch")
    n_sas <- sum(alpha_df$EthnicityTotal == "South-Asian Surinamese")

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
    write_csv(summary_table,
              paste0("results/alpha_diversity/alpha_diversity_summary_16s_",
                     site_name, "_ethnicity.csv"))

    ## ---- Wilcoxon rank-sum tests ----
    wilcox_results <- alpha_long |>
        group_by(metric) |>
        summarise(
            W        = wilcox.test(value ~ EthnicityTotal)$statistic,
            p.value  = wilcox.test(value ~ EthnicityTotal)$p.value,
            .groups  = "drop"
        )
    write_csv(wilcox_results,
              paste0("results/alpha_diversity/alpha_diversity_wilcoxon_16s_",
                     site_name, ".csv"))

    ## ---- Boxplots by ethnicity (primary figure) ----
    subtitle_text <- paste0("Dutch (n = ", n_dutch,
                            ") vs South-Asian Surinamese (n = ", n_sas, ")")

    ggplot(alpha_long, aes(x = EthnicityTotal, y = value, fill = EthnicityTotal)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
        stat_compare_means(method = "wilcox.test", label = "p.signif",
                           comparisons = list(c("Dutch", "South-Asian Surinamese"))) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = "Value", fill = "Ethnicity",
             title = paste0("Alpha diversity by ethnicity - 16S ", site_name),
             subtitle = subtitle_text,
             caption = paste0("Rarefied to ", format(rare_depth, big.mark = ","),
                              " reads/sample. Wilcoxon rank-sum test: ",
                              "ns p > 0.05, * p < 0.05, ** p < 0.01, ",
                              "*** p < 0.001, **** p < 0.0001")) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(paste0("results/alpha_diversity/alpha_diversity_boxplot_16s_",
                  site_name, "_ethnicity.pdf"),
           width = 10, height = 5)

    ## ---- Violin + boxplot (distribution overview) ----
    ggplot(alpha_long, aes(x = EthnicityTotal, y = value, fill = EthnicityTotal)) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.15, fill = "white",
                     outlier.shape = 21, outlier.size = 0.8) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = "Value", fill = "Ethnicity",
             title = paste0("Alpha diversity distribution - 16S ", site_name),
             subtitle = subtitle_text,
             caption = paste0("Rarefied to ", format(rare_depth, big.mark = ","),
                              " reads/sample")) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(paste0("results/alpha_diversity/alpha_diversity_violin_16s_",
                  site_name, "_ethnicity.pdf"),
           width = 10, height = 5)

    ## ---- Linear regression (adjusted for covariates) ----
    reg_results <- bind_rows(
        run_regression(alpha_df, "Observed", covariates),
        run_regression(alpha_df, "Shannon", covariates),
        run_regression(alpha_df, "Simpson", covariates)
    )
    write_csv(reg_results,
              paste0("results/alpha_diversity/alpha_diversity_regression_16s_",
                     site_name, ".csv"))

    cat("Completed:", site_name, "—", n_samples, "samples",
        "(Dutch:", n_dutch, "/ SAS:", n_sas, ")\n")
}
