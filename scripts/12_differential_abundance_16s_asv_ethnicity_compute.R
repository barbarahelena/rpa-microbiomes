## Compute ASV differential abundance: 16S throat/nose, pairwise ethnicity comparisons.
## Cohort selection, confounder screening and adjusted/unadjusted MaAsLin2 fits.
## Save model outputs and explicit reporting inputs; run the matching report
## to produce derived tables and plots. ASV and genus settings remain separate.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(pheatmap)

## Setup
setwd(here::here())

## Test mode: set DIFF_AB_TEST_N to cap each ethnicity group at N samples
## after group-size filtering, so the full pipeline runs in seconds instead
## of many minutes. Writes to a separate results dir so it can never clobber
## a real run. Example: DIFF_AB_TEST_N=40 Rscript scripts/12_differential_abundance_16s_asv_ethnicity_report.R
test_n <- suppressWarnings(as.integer(Sys.getenv("DIFF_AB_TEST_N", "")))
outdir <- if (!is.na(test_n)) "results/differential_abundance_test" else "results/differential_abundance"
if (!is.na(test_n)) cat("TEST MODE: capping each group at", test_n, "samples, writing to", outdir, "\n")

for (sub in c("cache", "maaslin2")) {
    dir.create(file.path(outdir, sub), recursive = TRUE, showWarnings = FALSE)
}

## Short abbreviations for filenames (group names contain spaces/hyphens)
eth_abbrev <- c(
    "Dutch" = "Dutch",
    "South-Asian Surinamese" = "SAS",
    "African Surinamese" = "AfrSur",
    "Javanese Surinamese" = "JavSur",
    "Other" = "Other",
    "Ghanaian" = "Ghanaian",
    "Turkish" = "Turkish",
    "Moroccan" = "Moroccan"
)
pair_name <- function(g1, g2) paste0(eth_abbrev[[g1]], "_vs_", eth_abbrev[[g2]])

## Keep only ethnicity groups with more than min_n samples at this site
keep_groups <- function(ps, min_n) {
    counts <- table(sample_data(ps)$EthnicityTotal)
    names(counts)[counts > min_n]
}

## Covariates to screen as potential confounders, per site. Screened once per
## site (see assess_confounders() below) across all qualifying ethnicity
## groups at once, so the resulting significant-confounder set is identical
## for every pairwise MaAsLin2 model at that site - not reassessed per pair.
## Site-specific lists follow the covariates that were significantly
## associated with beta diversity (Bray-Curtis and/or weighted UniFrac) for
## that site (see results/beta_diversity/covariate_screen_*_16s_<site>.csv).
## Antibiotics_FU is not screened here: antibiotic users are excluded from
## the analysis outright (see below) rather than adjusted for.
## MigrationGen and ResidenceDuration_BA are excluded: they are structurally
## NA for every Dutch participant, so they can never be part of a single
## covariate set that's shared across every pair - including the pairs that
## involve Dutch, which is exactly what "one consistent set per site"
## requires.
## DiscrMean_BA is excluded because it was only measured at baseline.
## Season is included at both sites: results/beta_diversity/covariate_screen/
## covariate_screen_*_16s_<site>.csv shows it significant (p=0.001) for both
## Bray-Curtis and weighted UniFrac, at both throat and nose - the largest R2
## of any covariate screened at nose, and second only to Smoking_FU at throat.
covariates_list <- list(
    throat = c("Age_FU", "Sex", "BMI_FU", "Smoking_FU", "AlcoholYN_FU",
               "HTSelfBP_FU", "DMSelfGluc_FU", "MetSyn_FU", "Lipidlowering_FU",
               "Antidepressants_FU", "Psychotropics_FU",
               "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU",
               "OralHealth_FU", "Season",
               "PM10_mean", "PM25_mean", "NO2_mean", "EC_mean"),
    nose = c("Age_FU", "Sex", "BMI_FU", "Smoking_FU", "AlcoholYN_FU",
             "MetSyn_FU", "Season",
             "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU",
             "PM10_mean", "PM25_mean", "NO2_mean", "EC_mean")
)

## Technical batch covariate, adjusted for in every MaAsLin2 model below
## regardless of significance (a technical artefact, not a candidate
## confounder to screen in/out) - so it's kept out of covariates_list/
## assess_confounders() and out of the "Adjusted for:" plot subtitles built
## from sig_confounders, same convention as always_covariates in
## 10_beta_diversity_16s_ethnicity_compute.R.
always_covariates <- c("SeqBatch")

## ---- Confounder assessment: omnibus across all qualifying groups ----
## Tests association between each covariate and ethnicity using every
## qualifying group at a site at once (not per pair), so a single set of
## significant confounders is shared across all of that site's pairwise
## MaAsLin2 models. Continuous covariates use Kruskal-Wallis (the >2-group
## generalization of Wilcoxon); categorical covariates use chi-square, or
## Fisher's exact (simulated p-value) when any expected cell count < 5 - the
## same contingency-table logic as before, already valid for any number of
## ethnicity groups.
assess_confounders <- function(meta, covariates) {
    lapply(covariates, function(cov) {
        vals <- meta[[cov]]
        eth <- meta$EthnicityTotal
        cc <- !is.na(vals)
        vals <- vals[cc]
        eth <- droplevels(factor(eth[cc]))

        ## Skip if too few complete cases, or if complete-case filtering
        ## leaves fewer than 2 ethnicity groups to compare
        if (sum(cc) < 10 || nlevels(eth) < 2) {
            return(tibble(covariate = cov, test = "skipped",
                          statistic = NA_real_, p.value = NA_real_))
        }

        if (is.factor(vals) || is.character(vals)) {
            vals <- droplevels(factor(vals))
            if (nlevels(vals) < 2) {
                return(tibble(covariate = cov, test = "skipped",
                              statistic = NA_real_, p.value = NA_real_))
            }
            tbl <- table(eth, vals)
            ## Use Fisher's exact test if any expected count < 5
            expected <- chisq.test(tbl)$expected
            if (any(expected < 5)) {
                res <- fisher.test(tbl, simulate.p.value = TRUE, B = 2000)
                return(tibble(covariate = cov, test = "fisher",
                              statistic = NA_real_, p.value = res$p.value))
            } else {
                res <- chisq.test(tbl)
                return(tibble(covariate = cov, test = "chisq",
                              statistic = res$statistic, p.value = res$p.value))
            }
        } else {
            res <- kruskal.test(vals ~ eth)
            return(tibble(covariate = cov, test = "kruskal",
                          statistic = res$statistic, p.value = res$p.value))
        }
    }) |> bind_rows()
}

## ---- Per-pair analysis: MaAsLin2 + visualization ----
## ps_site is already restricted to qualifying groups and antibiotic-free.
## sig_confounders was assessed once for the whole site (see
## assess_confounders() above) and is shared across every pair at this site.
## always_covariates (technical batch) is added to every model unconditionally
## and is never part of sig_confounders, so it never appears in the "Adjusted
## for:" plot subtitles below (which are built from sig_confounders only).
## Returns the existing summary row and the inputs needed by the report.
run_da_pair <- function(ps_site, site_name, group1, group2, sig_confounders,
                         always_covariates, outdir) {
    pair_tag <- pair_name(group1, group2)
    ## Short forms for plot titles only (axis/legend labels keep full names)
    g1_abbr <- eth_abbrev[[group1]]
    g2_abbr <- eth_abbrev[[group2]]
    cat("\n==", site_name, "-", group1, "vs", group2, "==\n")

    ## subset_samples() uses NSE that can't see group1/group2 when called
    ## from inside this function (its eval() looks one frame too high), so
    ## use prune_samples() with a plain logical vector instead
    keep_samples <- sample_names(ps_site)[
        as.character(sample_data(ps_site)$EthnicityTotal) %in% c(group1, group2)
    ]
    ps <- prune_samples(keep_samples, ps_site)

    ## Extract metadata; group1 is explicitly the reference level
    meta <- sample_data(ps) |>
        as("data.frame") |>
        mutate(EthnicityTotal = factor(EthnicityTotal, levels = c(group1, group2)))
    sample_data(ps) <- sample_data(meta)

    n_group1 <- sum(meta$EthnicityTotal == group1)
    n_group2 <- sum(meta$EthnicityTotal == group2)

    ## =========================================================================
    ## Step 2: MaAsLin2 differential abundance
    ## (sig_confounders is fixed for the whole site - see assess_confounders())
    ## =========================================================================

    ## Prepare input: taxa as columns, samples as rows
    counts_df <- as.data.frame(otu_table(ps))
    if (taxa_are_rows(ps)) {
        counts_df <- as.data.frame(t(counts_df))
    }

    ## Prepare metadata for MaAsLin2
    meta_maaslin <- meta |>
        select(EthnicityTotal, all_of(sig_confounders), all_of(always_covariates))

    ## Fixed effects: ethnicity + significant confounders + always_covariates
    ## (technical batch, adjusted for unconditionally - see always_covariates above)
    fixed_effects <- c("EthnicityTotal", sig_confounders, always_covariates)

    ## Output directory
    maaslin_outdir <- file.path(outdir, "maaslin2",
                                 paste0("maaslin2_16s_", site_name, "_", pair_tag))

    ## MaAsLin2's BH correction is computed independently within this model
    ## (across taxa). Unlike the ethnicity beta-diversity pairwise PERMANOVA - where every pair
    ## contributes exactly one p-value, making "all pairs" a well-defined BH
    ## family - each DA pair here produces its own family of per-taxon
    ## p-values, over a different sample subset and confounder set. Pooling
    ## q-values across pairs would conflate differently structured
    ## multiplicities, so FDR correction is deliberately kept independent per
    ## pair, matching how separate DESeq2/edgeR contrasts are each corrected
    ## on their own.
    maaslin_results <- Maaslin2(
        input_data     = counts_df,
        input_metadata = meta_maaslin,
        output         = maaslin_outdir,
        fixed_effects  = fixed_effects,
        normalization  = "TSS",
        transform      = "LOG",
        analysis_method = "LM",
        min_prevalence = 0.10,
        min_abundance  = 0.0001,
        reference      = paste0("EthnicityTotal,", group1),
        plot_heatmap   = FALSE,
        plot_scatter   = FALSE,
        max_significance = 0.25
    )

    ## Read results and filter to ethnicity effect
    res_all <- read_tsv(file.path(maaslin_outdir, "all_results.tsv"),
                        show_col_types = FALSE)
    res_eth <- res_all |>
        filter(metadata == "EthnicityTotal",
               value == group2)

    ## Get taxonomy lookup
    tax <- as.data.frame(tax_table(ps)) |>
        rownames_to_column("feature")
    res_eth <- res_eth |>
        left_join(tax, by = "feature")

    ## ---- Unadjusted model: ethnicity only, no confounders ----
    ## Fitted for comparison in the forest plot below.
    meta_maaslin_unadj <- meta |> select(EthnicityTotal)

    maaslin_outdir_unadj <- file.path(outdir, "maaslin2",
                                       paste0("maaslin2_16s_", site_name, "_",
                                              pair_tag, "_unadjusted"))

    maaslin_results_unadj <- Maaslin2(
        input_data     = counts_df,
        input_metadata = meta_maaslin_unadj,
        output         = maaslin_outdir_unadj,
        fixed_effects  = "EthnicityTotal",
        normalization  = "TSS",
        transform      = "LOG",
        analysis_method = "LM",
        min_prevalence = 0.10,
        min_abundance  = 0.0001,
        reference      = paste0("EthnicityTotal,", group1),
        plot_heatmap   = FALSE,
        plot_scatter   = FALSE,
        max_significance = 0.25
    )

    res_eth_unadj <- read_tsv(file.path(maaslin_outdir_unadj, "all_results.tsv"),
                              show_col_types = FALSE) |>
        filter(metadata == "EthnicityTotal", value == group2) |>
        select(feature, coef, stderr, pval, qval)

    n_sig <- sum(res_eth$qval < 0.25, na.rm = TRUE)
    n_strict <- sum(res_eth$qval < 0.05, na.rm = TRUE)
    cat(site_name, pair_tag, "- DA taxa (q < 0.25):", n_sig,
        "/ (q < 0.05):", n_strict, "\n")

    ## ggrepel formerly consumed this draw while rendering the volcano,
    ## between this pair's models and the next computation. Preserve that
    ## RNG sequence and pass the same draw to the separate report.
    volcano_seed <- if (n_sig > 0) sample.int(.Machine$integer.max, 1L) else NA_integer_

    cat("Completed:", site_name, pair_tag, "\n")

    list(
        summary = tibble(site = site_name, group1 = group1, group2 = group2,
            n_group1 = n_group1, n_group2 = n_group2,
            n_total = n_group1 + n_group2),
        report = list(
            ps = ps,
            meta = meta,
            site_name = site_name,
            group1 = group1,
            group2 = group2,
            sig_confounders = sig_confounders,
            pair_tag = pair_tag,
            g1_abbr = g1_abbr,
            g2_abbr = g2_abbr,
            n_group1 = n_group1,
            n_group2 = n_group2,
            res_eth = res_eth,
            res_eth_unadj = res_eth_unadj,
            n_sig = n_sig,
            n_strict = n_strict,
            volcano_seed = volcano_seed
        )
    )
}

## ---- Analysis loop over sites ----
sites <- list(
    throat = readRDS("data/processed/ps_throat_rarefied.RDS"),
    nose   = readRDS("data/processed/ps_nose_rarefied.RDS")
)

## Same N>50 threshold as the alpha/beta diversity analyses, so every
## ethnicity group screened there is also tested here.
min_group_n <- c(throat = 50, nose = 50)
summary_rows <- list()
report_pairs <- list()
confounders <- list()

for (site_name in names(sites)) {
    ps <- sites[[site_name]]
    covariates <- covariates_list[[site_name]]
    site_min_n <- min_group_n[[site_name]]

    ## Keep only ethnicity groups with N > site_min_n in this site
    qualifying <- keep_groups(ps, min_n = site_min_n)
    ps <- subset_samples(ps, EthnicityTotal %in% qualifying)

    ## Exclude participants on antibiotics (rather than adjusting for it)
    n_before_abx <- nsamples(ps)
    ps <- subset_samples(ps, Antibiotics_FU != "Yes")
    cat(site_name, "- excluded", n_before_abx - nsamples(ps),
        "antibiotic users\n")

    ## Extract metadata and drop unused factor levels
    meta <- sample_data(ps) |>
        as("data.frame") |>
        mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)))
    sample_data(ps) <- sample_data(meta)

    ## Test mode: cap each group at test_n samples (group eligibility above
    ## was already decided from the full data - every group here has >
    ## min_group_n samples, so test_n is always <= the group size)
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

    if (length(qualifying) < 2) {
        cat(site_name, "- qualifying groups (N >", site_min_n, "):",
            paste(qualifying, collapse = ", "),
            "| fewer than 2 groups qualify; skipping pairwise",
            "differential abundance for this site\n\n")
        next
    }

    pairs <- combn(qualifying, 2, simplify = FALSE)
    cat(site_name, "- qualifying groups (N >", site_min_n, "):",
        paste(qualifying, collapse = ", "), "| pairs to run:", length(pairs), "\n")

    ## ---- Confounder assessment: once per site, across all qualifying
    ## groups at once, so every pairwise MaAsLin2 model at this site shares
    ## the same adjustment set ----
    confounder_results <- assess_confounders(meta, covariates)
    confounders[[site_name]] <- confounder_results

    sig_confounders <- confounder_results |>
        filter(p.value < 0.05) |>
        pull(covariate)

    cat(site_name, "- Significant confounders (all groups):",
        paste(sig_confounders, collapse = ", "), "\n")

    for (pair in pairs) {
        pair_result <- run_da_pair(ps, site_name, pair[1], pair[2], sig_confounders,
                        always_covariates, outdir)
        summary_rows[[length(summary_rows) + 1]] <- pair_result$summary
        report_pairs[[length(report_pairs) + 1]] <- pair_result$report
    }

    cat("Completed:", site_name, "(", length(pairs), "pairs )\n\n")
}

dir.create(file.path(outdir, "cache"), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(pairs = report_pairs, confounders = confounders,
             summary = bind_rows(summary_rows)),
        file.path(outdir, "cache", "da_asv.rds"))
