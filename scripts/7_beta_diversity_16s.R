## Beta diversity analysis: 16S microbiome (throat and nose)
## Stratified by ethnicity (Dutch vs South-Asian Surinamese)
## with PERMANOVA, covariate screening, and betadisper

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
dir.create("results/beta_diversity", recursive = TRUE, showWarnings = FALSE)

## Define ethnicity colours
eth_colours <- c("Dutch" = "#1F78B4", "South-Asian Surinamese" = "#E31A1C")

## Covariates to screen
## Note: MigrationGen and ResidenceDuration_BA are excluded because they are
## structurally NA for all Dutch participants (migration-specific variables).
covariates <- c(
    # Demographics
    "Age_FU", "Sex", "BMI_FU",

    # Cardiometabolic risk factors
    "Smoking_FU", "AlcoholYN_FU", "SBP_FU", "DBP_FU",
    "HTSelfBP_FU", "DMSelfGluc_FU", "MetSyn_FU",

    # Medication
    "Antibiotics_FU", "Antihypertensiva_FU", "Lipidlowering_FU",
    "Corticosteroids_FU", "SystemicSteroids_FU", "Antihistamines_FU",
    "DecongAllerg_FU", "Antidepressants_FU", "Psychotropics_FU",

    # Mouth and nose variables
    "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU",
    "OralHealth_FU", "Nasal_FU"
)

## ---- Analysis loop over sites ----
sites <- list(
    throat = readRDS("data/processed/ps_throat_rarefied.RDS"),
    nose   = readRDS("data/processed/ps_nose_rarefied.RDS")
)

for (site_name in names(sites)) {
    ps <- sites[[site_name]]

    ## Filter to Dutch and South-Asian Surinamese only
    ps <- subset_samples(ps, EthnicityTotal %in% c("Dutch", "South-Asian Surinamese"))

    ## Extract metadata and drop unused factor levels
    meta <- sample_data(ps) |>
        as("data.frame") |>
        mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)))
    sample_data(ps) <- sample_data(meta)

    n_samples <- nsamples(ps)
    n_dutch <- sum(meta$EthnicityTotal == "Dutch")
    n_sas <- sum(meta$EthnicityTotal == "South-Asian Surinamese")
    subtitle_text <- paste0("Dutch (n = ", n_dutch,
                            ") vs South-Asian Surinamese (n = ", n_sas, ")")

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
            EthnicityTotal = meta$EthnicityTotal
        )

        ## PCoA plot coloured by ethnicity
        ggplot(ord_df, aes(x = PCo1, y = PCo2, colour = EthnicityTotal)) +
            geom_point(alpha = 0.5, size = 1) +
            stat_ellipse(level = 0.95, linewidth = 0.8) +
            scale_colour_manual(values = eth_colours) +
            labs(x = paste0("PCo1 (", var_explained[1], "%)"),
                 y = paste0("PCo2 (", var_explained[2], "%)"),
                 colour = "Ethnicity",
                 title = paste0("PCoA - ", dist_name, " - 16S ", site_name),
                 subtitle = subtitle_text) +
            theme_Publication() +
            theme(legend.position = "bottom")
        ggsave(paste0("results/beta_diversity/pcoa_", dist_label, "_16s_",
                      site_name, ".pdf"),
               width = 7, height = 6)

        ## ---- PERMANOVA: ethnicity only ----
        permanova_eth <- adonis2(
            dist_mat ~ EthnicityTotal,
            data = meta,
            permutations = 999
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

        ## ---- Full PERMANOVA: ethnicity + significant covariates ----
        if (length(sig_covariates) > 0) {
            formula_full <- as.formula(
                paste("dist_mat ~ EthnicityTotal +",
                      paste(sig_covariates, collapse = " + "))
            )
            ## Use complete cases for all variables in the model
            model_vars <- c("EthnicityTotal", sig_covariates)
            cc_idx <- complete.cases(meta[, model_vars])
            meta_cc <- meta[cc_idx, ] |>
                mutate(across(where(is.factor), droplevels))
            dist_cc <- as.dist(as.matrix(dist_mat)[cc_idx, cc_idx])

            permanova_full <- adonis2(
                as.formula(paste("dist_cc ~ EthnicityTotal +",
                                 paste(sig_covariates, collapse = " + "))),
                data = meta_cc,
                permutations = 999
            )
        } else {
            permanova_full <- permanova_eth
        }

        ## ---- Betadisper: test homogeneity of dispersions ----
        betadisp <- betadisper(dist_mat, meta$EthnicityTotal)
        betadisp_test <- permutest(betadisp, permutations = 999)

        ## Betadisper boxplot
        disp_df <- data.frame(
            Distance = betadisp$distances,
            EthnicityTotal = meta$EthnicityTotal
        )
        ggplot(disp_df, aes(x = EthnicityTotal, y = Distance,
                            fill = EthnicityTotal)) +
            geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
            scale_fill_manual(values = eth_colours) +
            labs(x = NULL, y = "Distance to centroid", fill = "Ethnicity",
                 title = paste0("Betadisper - ", dist_name, " - 16S ", site_name),
                 subtitle = paste0("Permutest p = ",
                                   format.pval(betadisp_test$tab[["Pr(>F)"]][1],
                                               digits = 3))) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1))
        ggsave(paste0("results/beta_diversity/betadisper_", dist_label, "_16s_",
                      site_name, ".pdf"),
               width = 5, height = 5)

        ## ---- Save results tables ----
        ## PERMANOVA ethnicity-only
        permanova_eth_df <- as.data.frame(permanova_eth) |>
            rownames_to_column("term") |>
            mutate(model = "ethnicity_only", .before = 1)

        ## PERMANOVA full model
        permanova_full_df <- as.data.frame(permanova_full) |>
            rownames_to_column("term") |>
            mutate(model = "adjusted", .before = 1)

        permanova_results <- bind_rows(permanova_eth_df, permanova_full_df)
        write_csv(permanova_results,
                  paste0("results/beta_diversity/permanova_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## Covariate screening
        write_csv(covariate_screen,
                  paste0("results/beta_diversity/covariate_screen_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## Betadisper
        betadisp_df <- tibble(
            F_stat  = betadisp_test$tab$F[1],
            p.value = betadisp_test$tab[["Pr(>F)"]][1]
        )
        write_csv(betadisp_df,
                  paste0("results/beta_diversity/betadisper_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## ---- Supplementary: PCoA coloured by significant covariates ----
        for (cov in sig_covariates) {
            ord_df[[cov]] <- meta[[cov]]
            p <- ggplot(ord_df, aes(x = PCo1, y = PCo2)) +
                geom_point(aes(colour = .data[[cov]]), alpha = 0.5, size = 1) +
                labs(x = paste0("PCo1 (", var_explained[1], "%)"),
                     y = paste0("PCo2 (", var_explained[2], "%)"),
                     colour = cov,
                     title = paste0("PCoA - ", dist_name, " - 16S ", site_name),
                     subtitle = paste0("Coloured by ", cov)) +
                theme_Publication() +
                theme(legend.position = "right")

            ## Use viridis for continuous, default for categorical
            if (is.numeric(meta[[cov]])) {
                p <- p + scale_colour_viridis_c(option = "plasma")
            }

            ggsave(paste0("results/beta_diversity/pcoa_", dist_label, "_16s_",
                          site_name, "_", cov, ".pdf"),
                   plot = p, width = 7, height = 6)
        }

        cat("Completed:", site_name, "-", dist_name, "\n")
    }

    cat("Finished site:", site_name, "-", n_samples, "samples\n\n")
}
