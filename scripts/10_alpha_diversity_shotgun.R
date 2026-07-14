## Alpha diversity analysis: shotgun metagenomics (tongue and throat)
## Stratified by ethnicity with linear regression

## Libraries
library(here)
library(tidyverse)
library(vegan)
library(broom)

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

## Define ethnicity colours
eth_colours <- c("Dutch" = "#1F78B4", "South-Asian Surinamese" = "#E31A1C")

## Covariates for linear regression
## Note: MigrationGen and ResidenceDuration_BA are excluded because they are
## structurally NA for all Dutch participants (migration-specific variables),
## which would drop the entire Dutch group from complete-case analysis.
covariates <- c("Age_FU", "Sex", "BMI_FU", "Smoking_FU", "Antibiotics_FU",
                "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU")

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
    write_csv(summary_table,
              paste0("results/alpha_diversity/alpha_diversity_summary_shotgun_",
                     site_name, ".csv"))

    ## ---- Boxplots by ethnicity (primary figure) ----
    ggplot(alpha_long, aes(x = EthnicityTotal, y = value, fill = EthnicityTotal)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = "Value", fill = "Ethnicity",
             title = paste0("Alpha diversity by ethnicity - shotgun ",
                            site_name, " (n = ", n_samples, ")")) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(paste0("results/alpha_diversity/alpha_diversity_boxplot_shotgun_",
                  site_name, ".pdf"),
           width = 8, height = 5)

    ## ---- Violin + boxplot (distribution overview) ----
    ggplot(alpha_long, aes(x = metric, y = value)) +
        geom_violin(fill = "#A6CEE3", alpha = 0.7) +
        geom_boxplot(width = 0.15, fill = "white",
                     outlier.shape = 21, outlier.size = 0.8) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        labs(x = NULL, y = "Value",
             title = paste0("Alpha diversity - shotgun ", site_name,
                            " (n = ", n_samples, ")")) +
        theme_Publication() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    ggsave(paste0("results/alpha_diversity/alpha_diversity_violin_shotgun_",
                  site_name, ".pdf"),
           width = 8, height = 5)

    ## ---- Linear regression ----
    reg_results <- bind_rows(
        run_regression(alpha_df, "Observed", covariates),
        run_regression(alpha_df, "Shannon", covariates)
    )
    write_csv(reg_results,
              paste0("results/alpha_diversity/alpha_diversity_regression_shotgun_",
                     site_name, ".csv"))

    cat("Completed:", site_name, "—", n_samples, "samples\n")
}
