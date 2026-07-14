## Beta diversity analysis: 16S microbiome (throat and nose)
## All non-Dutch ethnicity groups with N > 50, pooled (Dutch has no
## migration/acculturation data - see below)
## Screening migration/acculturation covariates (structurally NA for Dutch
## participants, so not usable in the ethnicity-comparison scripts)
## with PERMANOVA, covariate screening, and betadisper
## Primary grouping variable: MigrationGen (1st vs 2nd generation).
## Ethnicity is forced into every model as a covariate (entered before
## MigrationGen, so its variance is partialled out first): pooling several
## ethnic groups means their generation composition differs (e.g. Ghanaian
## is almost entirely 1st generation, n=2 2nd-gen), so without adjustment
## an apparent "generation" effect could really just be an ethnicity effect.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)

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
dir.create("results/beta_diversity_migration", recursive = TRUE, showWarnings = FALSE)

## Define migration generation colours
gen_colours <- c("1st generation" = "#33A02C", "2nd generation" = "#6A3D9A")

## Ethnicity colours (same validated palette as scripts 6-8, so a given
## ethnicity is always the same colour across every figure)
eth_colours <- c(
    "Dutch"                  = "#1F78B4",
    "South-Asian Surinamese" = "#E31A1C",
    "African Surinamese"     = "#33A02C",
    "Javanese Surinamese"    = "#6A3D9A",
    "Other"                  = "#B15928",
    "Ghanaian"               = "#FF7F00",
    "Turkish"                = "#E7298A",
    "Moroccan"               = "#D4AC0D"
)

## Keep only ethnicity groups with more than n=50 samples (matches Table 1)
keep_groups <- function(ps, min_n = 50) {
    counts <- table(sample_data(ps)$EthnicityTotal)
    names(counts)[counts > min_n]
}

## Migration/acculturation covariates to screen.
## These are structurally NA for all Dutch participants (baseline-only,
## migrant-specific variables), so they cannot be used in the Dutch vs
## South-Asian Surinamese comparison scripts. Here they are screened within
## every non-Dutch group with N > 50, pooled.
## Note: ResidenceDuration_BA and AgeMigration_BA are additionally NA for all
## 2nd-generation participants (born in the Netherlands, no migration event),
## so their screen is implicitly restricted to 1st-generation participants
## via complete-case filtering.
covariates <- c(
    "ResidenceDuration_BA", "AgeMigration_BA", "DifficultyDutch_BA",
    "CultFeelBerrys_BA", "CultOrientBerrys_BA", "CultNetworkBerrys_BA",
    "CultDistMeanScore0_BA", "CultDistMeanScore6_BA", "DiscrMean_BA"
)

## ---- Analysis loop over sites ----
sites <- list(
    throat = readRDS("data/processed/ps_throat_rarefied.RDS"),
    nose   = readRDS("data/processed/ps_nose_rarefied.RDS")
)

for (site_name in names(sites)) {
    ps <- sites[[site_name]]

    ## Filter to non-Dutch ethnicity groups with N > 50 (Dutch is excluded by
    ## construction anyway - MigrationGen is NA for every Dutch participant)
    qualifying <- setdiff(keep_groups(ps), "Dutch")
    ps <- subset_samples(ps, EthnicityTotal %in% qualifying)

    ## Extract metadata and drop unused factor levels
    meta <- sample_data(ps) |>
        as("data.frame") |>
        mutate(MigrationGen = droplevels(factor(MigrationGen)),
               EthnicityTotal = droplevels(factor(EthnicityTotal)))
    sample_data(ps) <- sample_data(meta)

    n_samples <- nsamples(ps)
    n_1st <- sum(meta$MigrationGen == "1st generation")
    n_2nd <- sum(meta$MigrationGen == "2nd generation")
    subtitle_text <- paste0("1st generation (n = ", n_1st,
                            ") vs 2nd generation (n = ", n_2nd, ")")
    gen_by_eth <- meta |> count(EthnicityTotal, MigrationGen)
    cat("Groups (N>50, non-Dutch) for", site_name, ":",
        paste(qualifying, collapse = ", "), "\n")
    cat("Generation split by ethnicity:\n")
    print(gen_by_eth)

    ## ---- Compute distance matrices ----
    dist_bc  <- phyloseq::distance(ps, method = "bray")
    dist_uni <- phyloseq::distance(ps, method = "wunifrac")

    distances <- list("Bray-Curtis" = dist_bc, "Weighted UniFrac" = dist_uni)

    for (dist_name in names(distances)) {
        dist_mat <- distances[[dist_name]]
        dist_label <- tolower(gsub("[- ]", "_", dist_name))

        ## ---- PCoA ordination ----
        pcoa <- ordinate(ps, method = "PCoA", distance = dist_mat)
        eig <- pcoa$values$Eigenvalues
        var_explained <- round(100 * eig / sum(eig), 1)

        ## Build ordination data frame
        ord_df <- data.frame(
            PCo1 = pcoa$vectors[, 1],
            PCo2 = pcoa$vectors[, 2],
            MigrationGen = meta$MigrationGen,
            EthnicityTotal = meta$EthnicityTotal
        )

        ## PCoA plot coloured by migration generation
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

        ## ---- PERMANOVA: ethnicity + migration generation ----
        ## Ethnicity enters first so MigrationGen's row reflects its effect
        ## net of ethnicity, not the ethnicity-confounded raw association.
        permanova_gen <- adonis2(
            dist_mat ~ EthnicityTotal + MigrationGen,
            data = meta,
            permutations = 999,
            by = "terms"
        )

        ## ---- Covariate screening (individual PERMANOVA per covariate) ----
        covariate_screen <- lapply(covariates, function(cov) {
            ## Use complete cases for this covariate
            cc_idx <- !is.na(meta[[cov]])
            if (sum(cc_idx) < 10) return(NULL)

            meta_cc <- meta[cc_idx, ] |>
                mutate(across(where(is.factor), droplevels))

            ## Check covariate has >= 2 levels
            vals <- meta_cc[[cov]]
            if (is.factor(vals) && nlevels(vals) < 2) return(NULL)
            if (!is.factor(vals) && length(unique(vals)) < 2) return(NULL)

            dist_cc <- as.dist(as.matrix(dist_mat)[cc_idx, cc_idx])

            formula <- as.formula(paste("dist_cc ~", cov))
            res <- adonis2(formula, data = meta_cc, permutations = 999)
            tibble(
                covariate = cov,
                Df        = res$Df[1],
                R2        = res$R2[1],
                F_stat    = res$F[1],
                p.value   = res[["Pr(>F)"]][1]
            )
        }) |> bind_rows()

        ## Identify significant covariates
        sig_covariates <- covariate_screen |>
            filter(p.value < 0.05) |>
            pull(covariate)

        ## ---- Full PERMANOVA: ethnicity + migration generation + significant
        ## covariates (ethnicity is forced, not screened - see header) ----
        ## ResidenceDuration_BA/AgeMigration_BA are structurally NA for every
        ## 2nd-generation participant (see header), so complete-case
        ## filtering on either together with MigrationGen drops all 2nd-gen
        ## rows, leaving MigrationGen with a single level - adonis2 can't fit
        ## a contrast for that. Keep them out of the joint model even if
        ## individually significant; their own univariate result (1st-gen
        ## only, by construction) is already saved in covariate_screen.
        model_covariates <- setdiff(sig_covariates,
                                    c("ResidenceDuration_BA", "AgeMigration_BA"))

        ## Any other covariate could in principle also happen to be NA for an
        ## entire MigrationGen or EthnicityTotal level in this particular
        ## complete-case subset - guard generically, not just for the two
        ## known offenders above, so a 45-min run never dies on this again.
        if (length(model_covariates) > 0) {
            model_vars <- c("EthnicityTotal", "MigrationGen", model_covariates)
            cc_idx <- complete.cases(meta[, model_vars])
            meta_cc <- meta[cc_idx, ] |>
                mutate(across(where(is.factor), droplevels))

            if (nlevels(meta_cc$EthnicityTotal) < 2 || nlevels(meta_cc$MigrationGen) < 2) {
                warning(site_name, " ", dist_name,
                        ": complete-case filtering on ", paste(model_covariates, collapse = ", "),
                        " collapsed EthnicityTotal or MigrationGen to a single level - ",
                        "falling back to the unadjusted model")
                permanova_full <- permanova_gen
            } else {
                dist_cc <- as.dist(as.matrix(dist_mat)[cc_idx, cc_idx])
                permanova_full <- adonis2(
                    as.formula(paste("dist_cc ~ EthnicityTotal + MigrationGen +",
                                     paste(model_covariates, collapse = " + "))),
                    data = meta_cc,
                    permutations = 999,
                    by = "terms"
                )
            }
        } else {
            permanova_full <- permanova_gen
        }

        ## ---- Betadisper: test homogeneity of dispersions ----
        betadisp <- betadisper(dist_mat, meta$MigrationGen)
        betadisp_test <- permutest(betadisp, permutations = 999)

        ## Betadisper boxplot
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

        cat("Completed:", site_name, "-", dist_name, "\n")
    }

    cat("Finished site:", site_name, "-", n_samples, "samples\n\n")
}
