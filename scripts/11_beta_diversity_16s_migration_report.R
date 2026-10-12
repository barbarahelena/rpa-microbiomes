## Report 16S beta diversity by migration generation (throat and nose).
## 1. Load the matching computation caches for pooled non-Dutch ethnicities N > 50.
## 2. Plot PCoA by generation, ethnicity, and significant covariates, plus dispersion.
## 3. Export cached PERMANOVA, covariate-screening, and dispersion tests;
##    PERMANOVA results account for ethnicity before migration generation.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)

## Functions
source(here::here("scripts", "lib", "plot_style.R"))

## Setup
setwd(here::here())
dir.create("results/beta_diversity_migration", recursive = TRUE, showWarnings = FALSE)

## Define migration generation colours
gen_colours <- c("1st generation" = "#33A02C", "2nd generation" = "#6A3D9A")

## Ethnicity colours (shared with alpha diversity, so a given
## ethnicity is always the same colour across every figure)
eth_colours <- ethnicity_colours

for (site_name in c("throat", "nose")) {
    for (dist_name in c("Bray-Curtis", "Weighted UniFrac")) {
        dist_label <- tolower(gsub("[- ]", "_", dist_name))
        ## ---- Load cached ordinations and tests for this site/distance ----
        cache_path <- file.path("results/beta_diversity_migration/cache", paste0("migration_", site_name, "_", dist_label, ".rds"))
        if (!file.exists(cache_path)) stop("Missing cache: ", cache_path, "; run pixi run beta-16s-migration-compute first.")
        cache <- readRDS(cache_path)
        meta <- cache$meta
        n_samples <- cache$n_samples
        subtitle_text <- cache$subtitle_text
        pcoa <- cache$pcoa
        var_explained <- cache$var_explained
        ord_df <- cache$ord_df
        permanova_gen <- cache$permanova_gen
        covariate_screen <- cache$covariate_screen
        sig_covariates <- cache$sig_covariates
        permanova_full <- cache$permanova_full
        betadisp <- cache$betadisp
        betadisp_test <- cache$betadisp_test
        ## ---- Plot PCoA by migration generation ----
        ggplot(ord_df, aes(x = PCo1, y = PCo2, colour = MigrationGen)) +
            geom_point(alpha = 0.5, size = 1) +
            stat_ellipse(level = 0.95, linewidth = 0.8) +
            scale_colour_manual(values = gen_colours) +
            labs(x = paste0("PCo1 (", var_explained[1], "%)"),
                 y = paste0("PCo2 (", var_explained[2], "%)"),
                 colour = "Migration generation",
                 title = paste0("PCoA - ", dist_name, " - 16S ", site_name,
                                " (all non-Dutch groups, N>50)"),
                 subtitle = subtitle_text) +
            theme_Publication() +
            theme(legend.position = "bottom")
        ggsave(paste0("results/beta_diversity_migration/pcoa_", dist_label, "_16s_",
                      site_name, ".pdf"),
               width = 7, height = 6)

        ## ---- Supplementary: same PCoA coloured by ethnicity, to see how
        ## much of the ordination it's driving before/next to MigrationGen ----
        ggplot(ord_df, aes(x = PCo1, y = PCo2, colour = EthnicityTotal)) +
            geom_point(alpha = 0.5, size = 1) +
            stat_ellipse(level = 0.95, linewidth = 0.8) +
            scale_colour_manual(values = eth_colours) +
            labs(x = paste0("PCo1 (", var_explained[1], "%)"),
                 y = paste0("PCo2 (", var_explained[2], "%)"),
                 colour = "Ethnicity",
                 title = paste0("PCoA - ", dist_name, " - 16S ", site_name,
                                " (all non-Dutch groups, N>50)")) +
            guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
            theme_Publication() +
            theme(legend.position = "bottom")
        ggsave(paste0("results/beta_diversity_migration/pcoa_", dist_label, "_16s_",
                      site_name, "_ethnicity.pdf"),
               width = 7, height = 8)

        ## ---- Plot dispersion by migration generation ----
        disp_df <- data.frame(
            Distance = betadisp$distances,
            MigrationGen = meta$MigrationGen
        )
        ggplot(disp_df, aes(x = MigrationGen, y = Distance,
                            fill = MigrationGen)) +
            geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
            scale_fill_manual(values = gen_colours) +
            labs(x = NULL, y = "Distance to centroid", fill = "Migration generation",
                 title = paste0("Betadisper - ", dist_name, " - 16S ", site_name,
                                " (all non-Dutch groups, N>50)"),
                 subtitle = paste0("Permutest p = ",
                                   format.pval(betadisp_test$tab[["Pr(>F)"]][1],
                                               digits = 3))) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1))
        ggsave(paste0("results/beta_diversity_migration/betadisper_", dist_label, "_16s_",
                      site_name, ".pdf"),
               width = 5, height = 5)

        ## ---- Save results tables ----
        ## PERMANOVA ethnicity + migration-generation (unadjusted for other covariates)
        permanova_gen_df <- as.data.frame(permanova_gen) |>
            rownames_to_column("term") |>
            mutate(model = "ethnicity_and_migrationgen", .before = 1)

        ## PERMANOVA full model
        permanova_full_df <- as.data.frame(permanova_full) |>
            rownames_to_column("term") |>
            mutate(model = "adjusted", .before = 1)

        permanova_results <- bind_rows(permanova_gen_df, permanova_full_df)
        write_csv(permanova_results,
                  paste0("results/beta_diversity_migration/permanova_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## Covariate screening
        write_csv(covariate_screen,
                  paste0("results/beta_diversity_migration/covariate_screen_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## Betadisper
        betadisp_df <- tibble(
            F_stat  = betadisp_test$tab$F[1],
            p.value = betadisp_test$tab[["Pr(>F)"]][1]
        )
        write_csv(betadisp_df,
                  paste0("results/beta_diversity_migration/betadisper_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## ---- Supplementary: PCoA coloured by significant covariates ----
        for (cov in sig_covariates) {
            ord_df[[cov]] <- meta[[cov]]
            p <- ggplot(ord_df, aes(x = PCo1, y = PCo2)) +
                geom_point(aes(colour = .data[[cov]]), alpha = 0.5, size = 1) +
                labs(x = paste0("PCo1 (", var_explained[1], "%)"),
                     y = paste0("PCo2 (", var_explained[2], "%)"),
                     colour = cov,
                     title = paste0("PCoA - ", dist_name, " - 16S ", site_name,
                                    " (all non-Dutch groups, N>50)"),
                     subtitle = paste0("Coloured by ", cov)) +
                theme_Publication() +
                theme(legend.position = "right")

            ## Use viridis for continuous, default for categorical
            if (is.numeric(meta[[cov]])) {
                p <- p + scale_colour_viridis_c(option = "plasma")
            }

            ggsave(paste0("results/beta_diversity_migration/pcoa_", dist_label, "_16s_",
                          site_name, "_", cov, ".pdf"),
                   plot = p, width = 7, height = 6)
        }

    }
}
