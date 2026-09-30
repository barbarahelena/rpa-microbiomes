## Report 16S alpha diversity by ethnicity (throat and nose; groups N > 50).
## 1. Load diversity values and test results from the matching computation cache.
## 2. Export summary statistics, group tests, and adjusted regression tables.
## 3. Save ethnicity boxplots with cached tests and violin plots by ethnicity.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(broom)
library(ggpubr)

## Functions
source(here::here("scripts", "lib", "plot_style.R"))

## Setup
setwd(here::here())
dir.create("results/alpha_diversity", recursive = TRUE, showWarnings = FALSE)

## Define ethnicity colours (shared with beta diversity, so a
## given ethnicity is always the same colour across every figure)
eth_colours <- ethnicity_colours

for (site_name in c("throat", "nose")) {
    ## ---- Load cached diversity values and results for this site ----
    cache_path <- paste0("results/alpha_diversity/cache/alpha_16s_", site_name, ".rds")
    if (!file.exists(cache_path)) stop("Missing cache: ", cache_path, "; run pixi run alpha-16s-compute first.")
    cache <- readRDS(cache_path)
    alpha_df <- cache$alpha_df
    alpha_long <- cache$alpha_long
    n_samples <- cache$n_samples
    summary_table <- cache$summary_table
    reg_results <- cache$reg_results
    group_ns <- cache$group_ns
    kruskal_results <- cache$kruskal_results
    pairwise_results <- cache$pairwise_results
    kruskal_annotations <- cache$kruskal_annotations
    figure_tests <- cache$figure_tests

    ## ---- Export summary statistics and test results ----
    write_csv(summary_table,
              paste0("results/alpha_diversity/alpha_diversity_summary_16s_",
                     site_name, "_ethnicity.csv"))
    write_csv(kruskal_results,
              paste0("results/alpha_diversity/alpha_diversity_kruskal_16s_",
                     site_name, ".csv"))
    write_csv(pairwise_results,
              paste0("results/alpha_diversity/alpha_diversity_pairwise_wilcoxon_16s_",
                     site_name, ".csv"))
    write_csv(reg_results,
              paste0("results/alpha_diversity/alpha_diversity_regression_16s_",
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
            step = (max_val - min_val) * 0.06,
            y.position = max_val + step * row_number(),
            p.adj.label = case_when(
                p.adj < 0.0001 ~ "****",
                p.adj < 0.001  ~ "***",
                p.adj < 0.01   ~ "**",
                TRUE           ~ "*"
            )
        ) |>
        ungroup()

    ## Original ggpubr label position: top of the trained panel scale,
    ## including bracket positions, with hjust = 0.2 and vjust = 0.
    annotation_positions <- metric_range |>
        left_join(sig_pairs |> group_by(metric) |>
                      summarise(bracket_top = max(y.position), .groups = "drop"),
                  by = "metric") |>
        mutate(y = pmax(max_val, bracket_top, na.rm = TRUE)) |>
        left_join(kruskal_annotations, by = "metric")

    p_box <- ggplot(alpha_long, aes(x = EthnicityTotal, y = value, fill = EthnicityTotal)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
        geom_text(data = annotation_positions,
                  aes(x = 1, y = y, label = paste0("p = ", p.format)),
                  inherit.aes = FALSE, hjust = 0.2, vjust = 0) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = "Value",
             title = paste0("Alpha diversity by ethnicity - 16S ", site_name)) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        theme_Publication() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 25, hjust = 1))

    if (nrow(sig_pairs) > 0) {
        p_box <- p_box +
            ggpubr::stat_pvalue_manual(sig_pairs, label = "p.adj.label",
                                        xmin = "group1", xmax = "group2",
                                        y.position = "y.position",
                                        tip.length = 0, size = 3,
                                        color = "grey30")
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
        labs(x = NULL, y = "Value",
             title = paste0("Alpha diversity distribution - 16S ", site_name)) +
        theme_Publication() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(paste0("results/alpha_diversity/alpha_diversity_violin_16s_",
                  site_name, "_ethnicity.pdf"),
           width = 12, height = 6.5)

}
