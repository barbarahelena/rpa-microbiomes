## Compute 16S beta diversity by migration generation (throat and nose).
## 1. Pool non-Dutch ethnicities with N > 50 from rarefied data; Dutch participants
##    lack migration/acculturation data.
## 2. Calculate Bray-Curtis/weighted UniFrac distances and PCoA ordinations.
## 3. Test migration generation, screen covariates, fit adjusted PERMANOVA,
##    and test dispersion; ethnicity enters first in every PERMANOVA model
##    to account for differing generation composition across ethnic groups.
## 4. Cache ordinations and test results for the matching report.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)

## Setup
setwd(here::here())
dir.create("results/beta_diversity_migration/cache", recursive = TRUE, showWarnings = FALSE)

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
## Note: ResidenceDuration_BA is additionally NA for all 2nd-generation
## participants (born in the Netherlands, no migration event), so its screen
## is implicitly restricted to 1st-generation participants via complete-case
## filtering. AgeMigration_BA is excluded here: it's more reflective of a
## participant's current age than of migration/acculturation. Only one of
## CultDistMeanScore0_BA/CultDistMeanScore6_BA is kept (they're the same
## underlying score on two different scales) to avoid a confusing duplicate.
covariates <- c(
    "ResidenceDuration_BA", "DifficultyDutch_BA",
    "CultFeelBerrys_BA", "CultOrientBerrys_BA", "CultNetworkBerrys_BA",
    "CultDistMeanScore0_BA", "DiscrMean_BA"
)

## ---- Load rarefied data and select non-Dutch ethnicity groups ----
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
        ## ResidenceDuration_BA is structurally NA for every 2nd-generation
        ## participant (see header), so complete-case filtering on it
        ## together with MigrationGen drops all 2nd-gen rows, leaving
        ## MigrationGen with a single level - adonis2 can't fit a contrast
        ## for that. Keep it out of the joint model even if individually
        ## significant; its own univariate result (1st-gen only, by
        ## construction) is already saved in covariate_screen.
        model_covariates <- setdiff(sig_covariates, "ResidenceDuration_BA")

        ## Any other covariate could in principle also happen to be NA for an
        ## entire MigrationGen or EthnicityTotal level in this particular
        ## complete-case subset - guard generically, not just for the known
        ## offender above, so a 45-min run never dies on this again.
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

        ## ---- Cache ordinations and test results ----
        saveRDS(list(meta = meta,
                n_samples = n_samples,
                subtitle_text = subtitle_text,
                pcoa = pcoa,
                var_explained = var_explained,
                ord_df = ord_df,
                permanova_gen = permanova_gen,
                covariate_screen = covariate_screen,
                sig_covariates = sig_covariates,
                permanova_full = permanova_full,
                betadisp = betadisp,
                betadisp_test = betadisp_test),
                file.path("results/beta_diversity_migration/cache",
                          paste0("migration_", site_name, "_", dist_label, ".rds")))

        cat("Completed:", site_name, "-", dist_name, "\n")
    }

    cat("Finished site:", site_name, "-", n_samples, "samples\n\n")
}
