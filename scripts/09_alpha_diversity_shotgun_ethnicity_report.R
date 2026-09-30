## Report shotgun alpha diversity by ethnicity (tongue and throat).
## 1. Load diversity values and results from the matching computation cache.
## 2. Export summary statistics and adjusted regression tables.
## 3. Save ethnicity boxplots and overall violin plots.

## Libraries
library(here)
library(tidyverse)
library(vegan)
library(broom)

## Functions
source(here::here("scripts", "lib", "plot_style.R"))

## Setup
setwd(here::here())
dir.create("results/alpha_diversity", recursive = TRUE, showWarnings = FALSE)

## Define ethnicity colours
eth_colours <- ethnicity_colours[c("Dutch", "South-Asian Surinamese")]

for (site_name in c("throat", "tongue")) {
    ## ---- Load cached diversity values and results for this site ----
    cache_path <- paste0("results/alpha_diversity/cache/alpha_shotgun_", site_name, ".rds")
    if (!file.exists(cache_path)) stop("Missing cache: ", cache_path, "; run pixi run alpha-shotgun-compute first.")
    cache <- readRDS(cache_path)
    alpha_df <- cache$alpha_df
    alpha_long <- cache$alpha_long
    n_samples <- cache$n_samples
    summary_table <- cache$summary_table
    reg_results <- cache$reg_results

    ## ---- Export summary statistics and test results ----
    write_csv(summary_table,
              paste0("results/alpha_diversity/alpha_diversity_summary_shotgun_",
                     site_name, ".csv"))
    write_csv(reg_results,
              paste0("results/alpha_diversity/alpha_diversity_regression_shotgun_",
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

}
