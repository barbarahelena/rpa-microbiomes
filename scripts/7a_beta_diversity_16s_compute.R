## Beta diversity analysis: 16S microbiome (throat and nose) - COMPUTE STEP
## Stratified by ethnicity (all groups with N > 50 per site)
## Runs every permutation-heavy PERMANOVA/betadisper call (the expensive part
## of this analysis) and caches the results to .rds so
## 7b_beta_diversity_16s_report.R can rebuild every plot/table without
## repeating the permutation tests.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)
library(parallel)

## Setup
setwd(here::here())

## Test mode: set BETA_DIV_TEST_N to cap each ethnicity group at N samples
## after group-size filtering, so the full pipeline runs in seconds instead
## of many minutes. Writes to a separate results dir so it can never clobber
## a real run. Example: BETA_DIV_TEST_N=40 Rscript scripts/7a_beta_diversity_16s_compute.R
test_n <- suppressWarnings(as.integer(Sys.getenv("BETA_DIV_TEST_N", "")))
outdir <- if (!is.na(test_n)) "results/beta_diversity_test" else "results/beta_diversity"
if (!is.na(test_n)) cat("TEST MODE: capping each group at", test_n, "samples, writing to", outdir, "\n")

dir.create(file.path(outdir, "cache"), recursive = TRUE, showWarnings = FALSE)

## Number of cores for parallelizing independent PERMANOVA calls (covariate
## screening, pairwise post-hoc, and the big omnibus/adjusted models via
## adonis2's own `parallel` arg). detectCores() reliably returns NA when
## sandboxed, so this is hardcoded for this 10-core machine (override with
## BETA_DIV_N_CORES) rather than relying on runtime detection.
n_cores <- suppressWarnings(as.integer(Sys.getenv("BETA_DIV_N_CORES", "8")))
if (is.na(n_cores)) n_cores <- 1

## Fail fast if the parallel setup is broken
self_test <- tryCatch({
    worker_out <- mclapply(1:4, function(i) i, mc.cores = n_cores)
    if (any(vapply(worker_out, function(x) inherits(x, "try-error"), logical(1))))
        stop("mclapply worker(s) failed")
    test_dist <- dist(matrix(rnorm(16), nrow = 4))
    test_grp  <- factor(c("a", "a", "b", "b"))
    adonis2(test_dist ~ test_grp, permutations = 9, parallel = n_cores)
    TRUE
}, error = function(e) e)
if (!isTRUE(self_test)) {
    stop("Parallel setup self-test failed (n_cores = ", n_cores, "): ",
         conditionMessage(self_test))
}
cat("Parallel setup OK - n_cores =", n_cores, "\n")

## Keep only ethnicity groups with more than n=50 samples (matches Table 1)
keep_groups <- function(ps, min_n = 50) { # this means categories such as Javanese Surinamese and Other are excluded
    counts <- table(sample_data(ps)$EthnicityTotal)
    names(counts)[counts > min_n]
}

## Pairwise PERMANOVA between every pair of groups (BH-adjusted).
## Pairs are independent, so run them across cores.
pairwise_permanova <- function(dist_mat, meta, group_var = "EthnicityTotal", n_cores = 1) {
    groups <- levels(droplevels(meta[[group_var]]))
    pairs <- combn(groups, 2, simplify = FALSE)
    mclapply(pairs, function(pair) {
        idx <- meta[[group_var]] %in% pair
        d_sub <- as.dist(as.matrix(dist_mat)[idx, idx])
        m_sub <- meta[idx, ] |> mutate(across(where(is.factor), droplevels))
        fit <- adonis2(as.formula(paste("d_sub ~", group_var)),
                        data = m_sub, permutations = 999)
        tibble(group1 = pair[1], group2 = pair[2],
               R2 = fit$R2[1], F_stat = fit$F[1], p.value = fit[["Pr(>F)"]][1])
    }, mc.cores = n_cores) |> bind_rows() |> mutate(p.adj = p.adjust(p.value, method = "BH"))
}

## Pairwise PERMANOVA between every pair of groups, adjusted for covariates
## (entered before the group term so the group's R2/p reflect its effect net
## of those covariates, matching the ethnicity-net-of-covariate convention
## used by the migration script's adjusted model). `covariates` is expected
## to be the significant covariates from this site x distance's screen (see
## sig_covariates in compute_permanova_block). Falls back to the unadjusted
## group-only model for a pair when complete-case filtering on the
## covariates collapses the group variable or a covariate to a single level
## within that pair.
pairwise_permanova_adjusted <- function(dist_mat, meta, covariates, unadjusted,
                                         group_var = "EthnicityTotal", n_cores = 1) {
    if (length(covariates) == 0) return(unadjusted)

    groups <- levels(droplevels(meta[[group_var]]))
    pairs <- combn(groups, 2, simplify = FALSE)
    mclapply(pairs, function(pair) {
        idx <- meta[[group_var]] %in% pair & complete.cases(meta[, covariates, drop = FALSE])
        m_sub <- meta[idx, ] |> mutate(across(where(is.factor), droplevels))

        group_ok <- nlevels(m_sub[[group_var]]) >= 2
        cov_ok <- vapply(covariates, function(cov) {
            v <- m_sub[[cov]]
            if (is.factor(v)) nlevels(v) >= 2 else length(unique(v)) >= 2
        }, logical(1))

        if (!group_ok || !all(cov_ok)) {
            idx_unadj <- meta[[group_var]] %in% pair
            d_unadj <- as.dist(as.matrix(dist_mat)[idx_unadj, idx_unadj])
            m_unadj <- meta[idx_unadj, ] |> mutate(across(where(is.factor), droplevels))
            fit <- adonis2(as.formula(paste("d_unadj ~", group_var)),
                            data = m_unadj, permutations = 999)
            return(tibble(group1 = pair[1], group2 = pair[2],
                           R2 = fit$R2[1], F_stat = fit$F[1],
                           p.value = fit[["Pr(>F)"]][1]))
        }

        d_sub <- as.dist(as.matrix(dist_mat)[idx, idx])
        fit <- adonis2(as.formula(paste("d_sub ~", paste(covariates, collapse = " + "),
                                        "+", group_var)),
                        data = m_sub, permutations = 999, by = "terms")
        tibble(group1 = pair[1], group2 = pair[2],
               R2 = fit[group_var, "R2"], F_stat = fit[group_var, "F"],
               p.value = fit[group_var, "Pr(>F)"])
    }, mc.cores = n_cores) |> bind_rows() |> mutate(p.adj = p.adjust(p.value, method = "BH"))
}

## Bundles every permutation-heavy PERMANOVA/betadisper call for one site x
## distance-metric combination into a single object.
compute_permanova_block <- function(dist_mat, meta, covariates, n_cores) {
    ## ---- PERMANOVA: ethnicity only (omnibus across all groups) ----
    permanova_eth <- adonis2(
        dist_mat ~ EthnicityTotal,
        data = meta,
        permutations = 999,
        parallel = n_cores,
        by = "terms"
    )

    ## ---- PERMANOVA: pairwise post-hoc between every pair of groups ----
    permanova_pairwise <- pairwise_permanova(dist_mat, meta, n_cores = n_cores)

    ## ---- Covariate screening (individual PERMANOVA per covariate) ----
    covariate_screen <- mclapply(covariates, function(cov) {
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
    }, mc.cores = n_cores) |> bind_rows()

    sig_covariates <- covariate_screen |>
        filter(p.value < 0.05) |>
        pull(covariate)

    ## ---- Ethnicity attenuation: how much does each covariate individually
    ## explain of ethnicity's effect on beta diversity? Fits
    ## dist ~ covariate + EthnicityTotal (covariate first, so ethnicity's R2
    ## is net of that one covariate) and compares it against ethnicity's R2 on
    ## that same covariate's complete-case subset (so the comparison isn't
    ## confounded by sample-size differences between covariates). A covariate
    ## that's a strong confounder of the ethnicity effect shows a large drop.
    ethnicity_attenuation <- mclapply(covariates, function(cov) {
        cc_idx <- !is.na(meta[[cov]])
        if (sum(cc_idx) < 10) return(NULL)

        meta_cc <- meta[cc_idx, ] |>
            mutate(across(where(is.factor), droplevels))

        vals <- meta_cc[[cov]]
        if (is.factor(vals) && nlevels(vals) < 2) return(NULL)
        if (!is.factor(vals) && length(unique(vals)) < 2) return(NULL)
        if (nlevels(meta_cc$EthnicityTotal) < 2) return(NULL)

        dist_cc <- as.dist(as.matrix(dist_mat)[cc_idx, cc_idx])

        eth_base <- adonis2(dist_cc ~ EthnicityTotal, data = meta_cc,
                             permutations = 999)
        fit <- adonis2(as.formula(paste("dist_cc ~", cov, "+ EthnicityTotal")),
                        data = meta_cc, permutations = 999, by = "terms")

        eth_R2_base <- eth_base$R2[1]
        eth_R2_net  <- fit["EthnicityTotal", "R2"]

        tibble(
            covariate               = cov,
            n                       = sum(cc_idx),
            eth_R2_unadjusted       = eth_R2_base,
            eth_R2_net_of_covariate = eth_R2_net,
            abs_reduction           = eth_R2_base - eth_R2_net,
            pct_reduction           = 100 * (eth_R2_base - eth_R2_net) / eth_R2_base,
            p_value                 = fit["EthnicityTotal", "Pr(>F)"]
        )
    }, mc.cores = n_cores) |> bind_rows() |> arrange(desc(abs_reduction))

    ## ---- Full PERMANOVA: ethnicity + significant covariates ----
    if (length(sig_covariates) > 0) {
        model_vars <- c("EthnicityTotal", sig_covariates)
        cc_idx <- complete.cases(meta[, model_vars])
        meta_cc <- meta[cc_idx, ] |>
            mutate(across(where(is.factor), droplevels))
        dist_cc <- as.dist(as.matrix(dist_mat)[cc_idx, cc_idx])

        permanova_full <- adonis2(
            as.formula(paste("dist_cc ~ EthnicityTotal +",
                             paste(sig_covariates, collapse = " + "))),
            data = meta_cc,
            permutations = 999,
            parallel = n_cores,
            by = "terms"
        )
    } else {
        permanova_full <- permanova_eth
    }

    ## ---- PERMANOVA pairwise post-hoc, adjusted for significant covariates ----
    permanova_pairwise_adjusted <- pairwise_permanova_adjusted(
        dist_mat, meta, sig_covariates, unadjusted = permanova_pairwise,
        n_cores = n_cores
    )

    ## ---- Betadisper: test homogeneity of dispersions (omnibus + pairwise) ----
    betadisp <- betadisper(dist_mat, meta$EthnicityTotal)
    betadisp_test <- permutest(betadisp, permutations = 999, parallel = n_cores)
    betadisp_pairwise <- as.data.frame(TukeyHSD(betadisp)$group) |>
        rownames_to_column("pair")

    list(
        permanova_eth               = permanova_eth,
        permanova_pairwise          = permanova_pairwise,
        permanova_pairwise_adjusted = permanova_pairwise_adjusted,
        covariate_screen            = covariate_screen,
        sig_covariates              = sig_covariates,
        ethnicity_attenuation       = ethnicity_attenuation,
        permanova_full              = permanova_full,
        betadisp                    = betadisp,
        betadisp_test               = betadisp_test,
        betadisp_pairwise           = betadisp_pairwise
    )
}

## Covariates to screen
## MigrationGen and ResidenceDuration_BA are excluded because NA for Dutch
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
    "OralHealth_FU", "Nasal_FU",

    # Air pollution exposure (RIVM ALO, 2013-2015 address-based mean)
    "PM10_mean", "PM25_mean", "NO2_mean", "EC_mean",

    # Sample collection
    "Season"
)

## ---- Compute loop over sites: throat and nose ----
sites <- list(
    throat = readRDS("data/processed/ps_throat_rarefied.RDS"),
    nose   = readRDS("data/processed/ps_nose_rarefied.RDS")
)

for (site_name in names(sites)) {
    ps <- sites[[site_name]]

    ## Filter to ethnicity groups with N > 50 in this site
    ps <- subset_samples(ps, EthnicityTotal %in% keep_groups(ps))

    ## Extract metadata and drop unused factor levels
    meta <- sample_data(ps) |>
        as("data.frame") |>
        mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)))
    sample_data(ps) <- sample_data(meta)

    ## Test mode: cap each group at test_n samples (group eligibility above
    ## was already decided from the full data - every group here has > 50
    ## samples, so test_n is always <= the group size)
    if (!is.na(test_n)) {
        set.seed(42)
        keep_samples <- meta |>
            rownames_to_column("sample_id") |>
            group_by(EthnicityTotal) |>
            slice_sample(n = test_n) |>
            pull(sample_id)
        ps <- prune_samples(keep_samples, ps)
        meta <- sample_data(ps) |>
            as("data.frame") |>
            mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)))
        sample_data(ps) <- sample_data(meta)
    }

    n_samples <- nsamples(ps)
    group_ns <- meta |> count(EthnicityTotal, name = "n") |> arrange(desc(n))
    cat("Groups (N>50) for", site_name, ":", paste(group_ns$EthnicityTotal, collapse = ", "), "\n")

    ## ---- Compute distance matrices ----
    distances <- list(
        "Bray-Curtis"      = phyloseq::distance(ps, method = "bray"),
        "Weighted UniFrac" = phyloseq::distance(ps, method = "wunifrac")
    )

    blocks <- list()
    pcoas <- list()

    for (dist_name in names(distances)) {
        dist_mat <- distances[[dist_name]]

        ## Every 999-permutation call for this site x distance combination is
        ## bundled into one object - the single most expensive step in this
        ## script.
        blocks[[dist_name]] <- compute_permanova_block(dist_mat, meta, covariates, n_cores)

        ## PCoA ordination - not permutation-based, but wunifrac requires the
        ## phylogenetic tree (ps), so it's cheapest to compute once here
        ## alongside the distance matrix rather than reloading ps downstream.
        pcoas[[dist_name]] <- ordinate(ps, method = "PCoA", distance = dist_mat)

        cat("Completed:", site_name, "-", dist_name, "\n")
    }

    saveRDS(
        list(
            meta      = meta,
            n_samples = n_samples,
            group_ns  = group_ns,
            distances = distances,
            pcoas     = pcoas,
            blocks    = blocks
        ),
        file.path(outdir, "cache", paste0("beta_diversity_16s_", site_name, ".rds"))
    )

    cat("Finished site:", site_name, "-", n_samples, "samples\n\n")
}
