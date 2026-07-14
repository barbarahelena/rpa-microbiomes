## Beta diversity analysis: 16S microbiome (throat and nose)
## Stratified by ethnicity (all groups with N > 50 per site)
## with PERMANOVA (omnibus + pairwise), covariate screening, and betadisper

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(vegan)
library(parallel)

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

## Test mode: set BETA_DIV_TEST_N to cap each ethnicity group at N samples
## after group-size filtering, so the full pipeline runs in seconds instead
## of many minutes. Writes to a separate results dir so it can never clobber
## a real run. Example: BETA_DIV_TEST_N=40 Rscript scripts/7_beta_diversity_16s.R
test_n <- suppressWarnings(as.integer(Sys.getenv("BETA_DIV_TEST_N", "")))
outdir <- if (!is.na(test_n)) "results/beta_diversity_test" else "results/beta_diversity"
if (!is.na(test_n)) cat("TEST MODE: capping each group at", test_n, "samples, writing to", outdir, "\n")

for (sub in c("pcoa", "permanova", "covariate_screen", "betadisper", "cache")) {
    dir.create(file.path(outdir, sub), recursive = TRUE, showWarnings = FALSE)
}

## Cache: the distance matrices and every 999-permutation PERMANOVA/betadisper
## call in this script are expensive, but most edits to this script only
## touch a downstream plot or table format. read_or_compute() loads a cached
## RDS instead of recomputing when one already exists at cache_file, so
## adding something new (like the pairwise heatmap below) doesn't require
## sitting through the full runtime again.
## IMPORTANT: this cache is NOT auto-invalidated if you change the group
## threshold, the covariates list, or the input data - delete <outdir>/cache/
## (or set BETA_DIV_FORCE_RECOMPUTE=1 for one run) after any such change.
force_recompute <- isTRUE(as.logical(Sys.getenv("BETA_DIV_FORCE_RECOMPUTE", "FALSE")))
if (force_recompute) cat("BETA_DIV_FORCE_RECOMPUTE=1 - ignoring any cached results\n")

## Cached results depend on test_n (a different sample subset each time it
## changes), so key the cache filename on it to avoid silently loading a
## stale cache computed under a different BETA_DIV_TEST_N value.
cache_suffix <- if (!is.na(test_n)) paste0("_testN", test_n) else ""

read_or_compute <- function(cache_file, compute_fn) {
    if (file.exists(cache_file) && !force_recompute) {
        cat("  [cache hit]", basename(cache_file), "\n")
        return(readRDS(cache_file))
    }
    result <- compute_fn()
    saveRDS(result, cache_file)
    result
}

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

## Define ethnicity colours
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

## Generic categorical palette for the supplementary covariate-coloured PCoA
## plots (same validated hues as eth_colours, applied by position since each
## covariate's factor levels differ)
cat_palette <- unname(eth_colours)

## Colours for the two distance metrics in the covariate-effect summary plot
dist_colours <- c("Bray-Curtis" = "#2a78d6", "Weighted UniFrac" = "#eb6834")

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

## Bundles every permutation-heavy PERMANOVA/betadisper call for one site x
## distance-metric combination into a single object, so it can be cached as
## one unit via read_or_compute() (see above).
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

    ## ---- Betadisper: test homogeneity of dispersions (omnibus + pairwise) ----
    betadisp <- betadisper(dist_mat, meta$EthnicityTotal)
    betadisp_test <- permutest(betadisp, permutations = 999, parallel = n_cores)
    betadisp_pairwise <- as.data.frame(TukeyHSD(betadisp)$group) |>
        rownames_to_column("pair")

    list(
        permanova_eth      = permanova_eth,
        permanova_pairwise = permanova_pairwise,
        covariate_screen   = covariate_screen,
        sig_covariates     = sig_covariates,
        permanova_full     = permanova_full,
        betadisp           = betadisp,
        betadisp_test      = betadisp_test,
        betadisp_pairwise  = betadisp_pairwise
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
    "OralHealth_FU", "Nasal_FU"
)

## Human-readable labels for plots (raw variable names stay in filenames/CSVs)
covariate_labels <- c(
    Age_FU               = "Age",
    Sex                  = "Sex",
    BMI_FU               = "BMI",
    Smoking_FU           = "Smoking status",
    AlcoholYN_FU         = "Alcohol use",
    SBP_FU               = "Systolic blood pressure",
    DBP_FU               = "Diastolic blood pressure",
    HTSelfBP_FU          = "Hypertension",
    DMSelfGluc_FU        = "Diabetes",
    MetSyn_FU            = "Metabolic syndrome",
    Antibiotics_FU       = "Antibiotics use",
    Antihypertensiva_FU  = "Blood pressure lowering drugs",
    Lipidlowering_FU     = "Lipid lowering drugs",
    Corticosteroids_FU   = "Corticosteroids use",
    SystemicSteroids_FU  = "Systemic steroids use",
    Antihistamines_FU    = "Antihistamines use",
    DecongAllerg_FU      = "Decongestant/allergy medication",
    Antidepressants_FU   = "Antidepressants use",
    Psychotropics_FU     = "Psychotropic medication use",
    ToothBrushing_FU     = "Tooth brushing frequency",
    TongueBrushing_FU    = "Tongue brushing frequency",
    Mouthwash_FU         = "Mouthwash use",
    OralHealth_FU        = "Self-rated oral health",
    Nasal_FU             = "Nasal medication use",
    EthnicityTotal       = "Ethnicity"
)

## ---- Analysis loop over sites: throat and nose ----
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
    distances <- read_or_compute(
        file.path(outdir, "cache", paste0("distances_16s_", site_name, cache_suffix, ".rds")),
        function() list(
            "Bray-Curtis"      = phyloseq::distance(ps, method = "bray"),
            "Weighted UniFrac" = phyloseq::distance(ps, method = "wunifrac")
        )
    )

    ## Collects ethnicity + covariate PERMANOVA R2/p from both distance
    ## metrics, for the covariate-effect summary plot built after this loop.
    covariate_screen_all <- list()

    for (dist_name in names(distances)) {
        dist_mat <- distances[[dist_name]]
        dist_label <- tolower(gsub("[- ]", "_", dist_name))

        ## ---- PERMANOVA (omnibus + pairwise), covariate screen, betadisper ----
        ## Every 999-permutation call for this site x distance combination is
        ## bundled into one cached object (see read_or_compute() above) -
        ## the single most expensive step in this script.
        block <- read_or_compute(
            file.path(outdir, "cache",
                      paste0("permanova_block_", dist_label, "_16s_",
                             site_name, cache_suffix, ".rds")),
            function() compute_permanova_block(dist_mat, meta, covariates, n_cores)
        )
        permanova_eth      <- block$permanova_eth
        permanova_pairwise <- block$permanova_pairwise
        covariate_screen   <- block$covariate_screen
        sig_covariates     <- block$sig_covariates
        permanova_full     <- block$permanova_full
        betadisp           <- block$betadisp
        betadisp_test      <- block$betadisp_test
        betadisp_pairwise  <- block$betadisp_pairwise

        permanova_label <- paste0(
            "PERMANOVA: R² = ", round(permanova_eth$R2[1], 3),
            ", p = ", format.pval(permanova_eth[["Pr(>F)"]][1], digits = 2, eps = 0.001)
        )

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
        ## Always show large 95% ellipses. 4+ groups additionally get small
        ## centroids on top (ellipses alone get hard to pin down once there
        ## are several overlapping groups), plus the omnibus PERMANOVA R2/p
        ## annotated top-right.
        n_groups <- nlevels(droplevels(meta$EthnicityTotal))

        p_pcoa <- ggplot(ord_df, aes(x = PCo1, y = PCo2, colour = EthnicityTotal)) +
            geom_point(alpha = 0.5, size = 1) +
            stat_ellipse(level = 0.95, linewidth = 0.8)

        if (n_groups > 3) {
            centroids <- ord_df |>
                group_by(EthnicityTotal) |>
                summarise(PCo1 = mean(PCo1), PCo2 = mean(PCo2), .groups = "drop")
            p_pcoa <- p_pcoa +
                geom_point(data = centroids,
                           aes(x = PCo1, y = PCo2, fill = EthnicityTotal),
                           shape = 21, colour = "black", size = 4, stroke = 0.8) +
                scale_fill_manual(values = eth_colours, guide = "none")
        }

        p_pcoa <- p_pcoa +
            scale_colour_manual(values = eth_colours) +
            annotate("text", x = Inf, y = Inf, label = permanova_label,
                     hjust = 1.05, vjust = 1.5, size = 3.2) +
            labs(x = paste0("PCo1 (", var_explained[1], "%)"),
                 y = paste0("PCo2 (", var_explained[2], "%)"),
                 colour = "Ethnicity",
                 title = paste0("PCoA - ", dist_name, " - 16S ", site_name)) +
            guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
            theme_Publication() +
            theme(legend.position = "bottom")
        ggsave(paste0(outdir, "/pcoa/pcoa_", dist_label, "_16s_",
                      site_name, ".pdf"),
               plot = p_pcoa, width = 7, height = 8)

        ## Stash ethnicity + covariate effect sizes for the summary plot
        covariate_screen_all[[dist_name]] <- bind_rows(
            tibble(covariate = "EthnicityTotal", Df = permanova_eth$Df[1],
                   R2 = permanova_eth$R2[1], F_stat = permanova_eth$F[1],
                   p.value = permanova_eth[["Pr(>F)"]][1]),
            covariate_screen
        ) |> mutate(distance = dist_name)

        ## Betadisper boxplot
        disp_df <- data.frame(
            Distance = betadisp$distances,
            EthnicityTotal = meta$EthnicityTotal
        )
        ggplot(disp_df, aes(x = EthnicityTotal, y = Distance,
                            fill = EthnicityTotal)) +
            geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
            scale_fill_manual(values = eth_colours, guide = "none") +
            labs(x = NULL, y = "Distance to centroid",
                 title = paste0("Betadisper - ", dist_name, " - 16S ", site_name),
                 subtitle = paste0("Permutest p = ",
                                   format.pval(betadisp_test$tab[["Pr(>F)"]][1],
                                               digits = 3))) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1))
        ggsave(paste0(outdir, "/betadisper/betadisper_", dist_label, "_16s_",
                      site_name, ".pdf"),
               width = 5, height = 5)

        ## ---- Save results tables ----
        ## PERMANOVA ethnicity-only and full (adjusted) model, joined by term
        ## so R2/p for the same term are side by side across both models.
        permanova_eth_df <- as.data.frame(permanova_eth) |>
            rownames_to_column("term") |>
            select(term, Df, R2, `F`, p.value = `Pr(>F)`) |>
            rename_with(~ paste0(., "_ethnicity_only"), -term)

        permanova_full_df <- as.data.frame(permanova_full) |>
            rownames_to_column("term") |>
            select(term, Df, R2, `F`, p.value = `Pr(>F)`) |>
            rename_with(~ paste0(., "_adjusted"), -term)

        permanova_results <- full_join(permanova_eth_df, permanova_full_df, by = "term")
        write_csv(permanova_results,
                  paste0(outdir, "/permanova/permanova_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## PERMANOVA pairwise post-hoc
        write_csv(permanova_pairwise,
                  paste0(outdir, "/permanova/permanova_pairwise_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## ---- Pairwise PERMANOVA heatmap (R2, BH-adjusted significance) ----
        ## Mirror every pair onto both triangles so geom_tile draws a full
        ## symmetric group x group grid instead of a half-empty one.
        groups_order <- levels(droplevels(meta$EthnicityTotal))
        pairwise_mat_df <- bind_rows(
            permanova_pairwise |> select(group1, group2, R2, p.adj),
            permanova_pairwise |> select(group1 = group2, group2 = group1, R2, p.adj)
        ) |>
            mutate(
                group1 = factor(group1, levels = groups_order),
                group2 = factor(group2, levels = groups_order),
                sig = case_when(
                    p.adj < 0.0001 ~ "****",
                    p.adj < 0.001  ~ "***",
                    p.adj < 0.01   ~ "**",
                    p.adj < 0.05   ~ "*",
                    TRUE           ~ ""
                ),
                cell_label = paste0(sprintf("%.3f", R2), "\n", sig)
            )

        p_permanova_heat <- ggplot(pairwise_mat_df, aes(x = group1, y = group2, fill = R2)) +
            geom_tile(colour = "white") +
            geom_text(aes(label = cell_label), size = 3, lineheight = 0.9) +
            scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", limits = c(0, NA)) +
            scale_x_discrete(drop = FALSE) +
            scale_y_discrete(drop = FALSE) +
            labs(x = NULL, y = NULL, fill = expression(R^2),
                 title = paste0("Pairwise PERMANOVA - ", dist_name, " - 16S ", site_name),
                 caption = "BH-adjusted: * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001") +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "right",
                  panel.grid = element_blank())
        ggsave(paste0(outdir, "/permanova/permanova_pairwise_heatmap_", dist_label,
                      "_16s_", site_name, ".pdf"),
               plot = p_permanova_heat, width = 6, height = 5.5)

        ## Covariate screening
        write_csv(covariate_screen,
                  paste0(outdir, "/covariate_screen/covariate_screen_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## Betadisper
        betadisp_df <- tibble(
            F_stat  = betadisp_test$tab$F[1],
            p.value = betadisp_test$tab[["Pr(>F)"]][1]
        )
        write_csv(betadisp_df,
                  paste0(outdir, "/betadisper/betadisper_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## Betadisper pairwise post-hoc (Tukey HSD)
        write_csv(betadisp_pairwise,
                  paste0(outdir, "/betadisper/betadisper_pairwise_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## ---- Supplementary: PCoA coloured by significant covariates ----
        for (cov in sig_covariates) {
            ord_df[[cov]] <- meta[[cov]]
            ## Drop cases with a missing value for this covariate before plotting
            plot_df <- ord_df[!is.na(ord_df[[cov]]), ]

            cov_stats <- covariate_screen |> filter(covariate == cov)
            cov_label <- paste0(
                "PERMANOVA: R² = ", round(cov_stats$R2[1], 3),
                ", p = ", format.pval(cov_stats$p.value[1], digits = 2, eps = 0.001)
            )

            p <- ggplot(plot_df, aes(x = PCo1, y = PCo2, colour = .data[[cov]])) +
                geom_point(alpha = 0.5, size = 1) +
                annotate("text", x = Inf, y = Inf, label = cov_label,
                         hjust = 1.05, vjust = 1.5, size = 3.2) +
                labs(x = paste0("PCo1 (", var_explained[1], "%)"),
                     y = paste0("PCo2 (", var_explained[2], "%)"),
                     colour = covariate_labels[[cov]],
                     title = paste0("PCoA - ", dist_name, " - 16S ", site_name)) +
                theme_Publication() +
                theme(legend.position = "right")

            if (is.numeric(plot_df[[cov]])) {
                ## Continuous covariate: viridis, no ellipse/centroid.
                ## theme_Publication()'s legend.key.size is tiny (sized for
                ## discrete dot legends), so give the colourbar its own height.
                p <- p + scale_colour_viridis_c(
                    option = "plasma",
                    guide = guide_colourbar(barheight = unit(4, "cm"))
                )
            } else {
                ## Categorical covariate: validated palette, always a large
                ## ellipse, plus centroids on top once there are 4+ levels -
                ## same rule as the main ethnicity plot
                p <- p +
                    scale_colour_manual(values = cat_palette) +
                    stat_ellipse(level = 0.95, linewidth = 0.8)
                n_cov_levels <- nlevels(droplevels(factor(plot_df[[cov]])))
                if (n_cov_levels > 3) {
                    cov_centroids <- plot_df |>
                        group_by(.data[[cov]]) |>
                        summarise(PCo1 = mean(PCo1), PCo2 = mean(PCo2), .groups = "drop")
                    p <- p +
                        geom_point(data = cov_centroids,
                                   aes(x = PCo1, y = PCo2, fill = .data[[cov]]),
                                   shape = 21, colour = "black", size = 4, stroke = 0.8) +
                        scale_fill_manual(values = cat_palette, guide = "none")
                }
            }

            ggsave(paste0(outdir, "/pcoa/pcoa_", dist_label, "_16s_",
                          site_name, "_", cov, ".pdf"),
                   plot = p, width = 7, height = 6)
        }

        cat("Completed:", site_name, "-", dist_name, "\n")
    }

    ## ---- Summary plot: ethnicity + covariate R2 across both distance metrics ----
    covariate_effects <- bind_rows(covariate_screen_all)

    ## Order covariates by effect size (ascending, so largest ends up at the
    ## top after coord_flip); EthnicityTotal always sits at the very top as
    ## the primary comparison, not ranked in among the covariates.
    covariate_order <- covariate_effects |>
        filter(covariate != "EthnicityTotal") |>
        group_by(covariate) |>
        summarise(max_R2 = max(R2), .groups = "drop") |>
        arrange(max_R2) |>
        pull(covariate)
    covariate_effects <- covariate_effects |>
        mutate(covariate_label = covariate_labels[covariate],
               covariate_label = factor(covariate_label,
                                        levels = covariate_labels[c(covariate_order, "EthnicityTotal")]),
               significant = p.value < 0.05)

    write_csv(covariate_effects,
              paste0(outdir, "/permanova/permanova_summary_16s_", site_name, ".csv"))

    ## Significant bars are filled solid; non-significant bars show only
    ## their coloured outline (fill alpha = 0) - no text label needed, and
    ## it sidesteps ever having to align a star with a dodged bar.
    n_terms <- n_distinct(covariate_effects$covariate_label)
    ggplot(covariate_effects, aes(x = covariate_label, y = R2, fill = distance,
                                   colour = distance, alpha = significant)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.6) +
        coord_flip() +
        scale_fill_manual(values = dist_colours) +
        scale_colour_manual(values = dist_colours, guide = "none") +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
        labs(x = NULL, y = expression(R^2), fill = "Distance metric",
             title = paste0("Covariate effects on beta diversity - 16S ", site_name),
             caption = "Filled bars: p < 0.05 (PERMANOVA)") +
        theme_Publication()
    ggsave(paste0(outdir, "/permanova/permanova_summary_16s_", site_name, ".pdf"),
           width = 8, height = max(6, 0.3 * n_terms + 2))

    cat("Finished site:", site_name, "-", n_samples, "samples\n\n")
}
