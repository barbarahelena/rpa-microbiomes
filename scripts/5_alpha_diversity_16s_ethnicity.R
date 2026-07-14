## Alpha diversity analysis: 16S microbiome (throat and nose)
## Stratified by ethnicity (all groups with N > 50 per site)
## with Kruskal-Wallis (+ pairwise Wilcoxon) and covariate-adjusted regression

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

## Define ethnicity colours (same validated palette as scripts 7 and 8, so a
## given ethnicity is always the same colour across every figure)
eth_colours <- c(
    "Dutch"                  = "#1F78B4",  # blue
    "South-Asian Surinamese" = "#E31A1C",  # red
    "African Surinamese"     = "#33A02C",  # green
    "Javanese Surinamese"    = "#6A3D9A",  # purple
    "Other"                  = "#B15928",  # brown
    "Ghanaian"               = "#FF7F00",  # orange
    "Turkish"                = "#E7298A",  # magenta
    "Moroccan"               = "#D4AC0D"   # gold
)

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
    write_csv(summary_table,
              paste0("results/alpha_diversity/alpha_diversity_summary_16s_",
                     site_name, "_ethnicity.csv"))

    ## ---- Kruskal-Wallis omnibus test across all groups ----
    kruskal_results <- alpha_long |>
        group_by(metric) |>
        summarise(
            statistic = kruskal.test(value ~ EthnicityTotal)$statistic,
            p.value   = kruskal.test(value ~ EthnicityTotal)$p.value,
            .groups   = "drop"
        )
    write_csv(kruskal_results,
              paste0("results/alpha_diversity/alpha_diversity_kruskal_16s_",
                     site_name, ".csv"))

    ## ---- Pairwise Wilcoxon rank-sum tests between every pair of groups ----
    pairwise_results <- alpha_long |>
        group_by(metric) |>
        group_modify(~ {
            pw <- pairwise.wilcox.test(.x$value, .x$EthnicityTotal, p.adjust.method = "BH")
            as.data.frame(as.table(pw$p.value)) |>
                filter(!is.na(Freq)) |>
                rename(group1 = Var1, group2 = Var2, p.adj = Freq)
        }) |>
        ungroup()
    write_csv(pairwise_results,
              paste0("results/alpha_diversity/alpha_diversity_pairwise_wilcoxon_16s_",
                     site_name, ".csv"))

    ## ---- Boxplots by ethnicity (primary figure) ----
    ## Kruskal-Wallis omnibus p-value per facet, plus brackets for pairwise
    ## Wilcoxon (BH-adjusted) comparisons with p.adj < 0.05 only - with up to
    ## 15 (throat) or 10 (nose) possible pairs, showing every pair would be
    ## unreadable. Full pairwise results (significant or not) are always in
    ## the CSV regardless of what gets plotted.
    metric_range <- alpha_long |>
        group_by(metric) |>
        summarise(max_val = max(value), min_val = min(value), .groups = "drop")

    sig_pairs <- pairwise_results |>
        filter(p.adj < 0.05) |>
        left_join(metric_range, by = "metric") |>
        group_by(metric) |>
        arrange(p.adj) |>
        mutate(
            group1 = as.character(group1),
            group2 = as.character(group2),
            step = (max_val - min_val) * 0.12,
            y.position = max_val + step * row_number(),
            p.adj.label = case_when(
                p.adj < 0.0001 ~ "****",
                p.adj < 0.001  ~ "***",
                p.adj < 0.01   ~ "**",
                TRUE           ~ "*"
            )
        ) |>
        ungroup()

    p_box <- ggplot(alpha_long, aes(x = EthnicityTotal, y = value, fill = EthnicityTotal)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
        stat_compare_means(method = "kruskal.test", label = "p.format") +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = "Value", fill = "Ethnicity",
             title = paste0("Alpha diversity by ethnicity - 16S ", site_name),
             caption = paste0("Rarefied to ", format(rare_depth, big.mark = ","),
                              " reads/sample. Kruskal-Wallis omnibus test across all groups; ",
                              "brackets show pairwise Wilcoxon (BH-adjusted) p < 0.05 only - ",
                              "full pairwise results in results CSV.")) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))

    if (nrow(sig_pairs) > 0) {
        p_box <- p_box +
            ggpubr::stat_pvalue_manual(sig_pairs, label = "p.adj.label",
                                        xmin = "group1", xmax = "group2",
                                        y.position = "y.position",
                                        tip.length = 0.01, size = 3)
    }

    ggsave(paste0("results/alpha_diversity/alpha_diversity_boxplot_16s_",
                  site_name, "_ethnicity.pdf"),
           plot = p_box, width = 12, height = 6.5)

    ## ---- Violin + boxplot (distribution overview) ----
    ggplot(alpha_long, aes(x = EthnicityTotal, y = value, fill = EthnicityTotal)) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.15, fill = "white",
                     outlier.shape = 21, outlier.size = 0.8) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = "Value", fill = "Ethnicity",
             title = paste0("Alpha diversity distribution - 16S ", site_name),
             caption = paste0("Rarefied to ", format(rare_depth, big.mark = ","),
                              " reads/sample")) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(paste0("results/alpha_diversity/alpha_diversity_violin_16s_",
                  site_name, "_ethnicity.pdf"),
           width = 12, height = 6.5)

    ## ---- Linear regression (adjusted for covariates) ----
    reg_results <- bind_rows(
        run_regression(alpha_df, "Observed", covariates),
        run_regression(alpha_df, "Shannon", covariates),
        run_regression(alpha_df, "Simpson", covariates)
    )
    write_csv(reg_results,
              paste0("results/alpha_diversity/alpha_diversity_regression_16s_",
                     site_name, ".csv"))

    cat("Completed:", site_name, "—", n_samples, "samples across",
        nrow(group_ns), "groups\n")
}
