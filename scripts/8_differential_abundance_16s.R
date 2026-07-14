## Differential abundance analysis: 16S microbiome (throat and nose)
## Pairwise ethnicity comparisons (all groups with N above a per-site
## minimum - throat requires N > 400, nose N > 150 since it has far fewer
## samples per group)
## Step 1: Confounder assessment (once per site, across all qualifying
##         groups at once, so every pairwise model at a site shares the
##         same adjustment set)
## Step 2: MaAsLin2 adjusted for significant confounders
## Step 3: Visualization (volcano, heatmap, taxon boxplots)

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(pheatmap)

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

## Test mode: set DIFF_AB_TEST_N to cap each ethnicity group at N samples
## after group-size filtering, so the full pipeline runs in seconds instead
## of many minutes. Writes to a separate results dir so it can never clobber
## a real run. Example: DIFF_AB_TEST_N=40 Rscript scripts/8_differential_abundance_16s.R
test_n <- suppressWarnings(as.integer(Sys.getenv("DIFF_AB_TEST_N", "")))
outdir <- if (!is.na(test_n)) "results/differential_abundance_test" else "results/differential_abundance"
if (!is.na(test_n)) cat("TEST MODE: capping each group at", test_n, "samples, writing to", outdir, "\n")

for (sub in c("confounder_assessment", "maaslin2", "volcano", "heatmap", "forestplot", "boxplots")) {
    dir.create(file.path(outdir, sub), recursive = TRUE, showWarnings = FALSE)
}

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

## Shorten covariate names for plot subtitles (drop trailing _FU/_BA suffix)
## and wrap long lists onto multiple lines instead of running off the page.
format_covariates <- function(covs, width = 70) {
    covs |>
        gsub("_(FU|BA)$", "", x = _) |>
        paste(collapse = ", ") |>
        str_wrap(width = width)
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
covariates_list <- list(
    throat = c("Age_FU", "Sex", "BMI_FU", "Smoking_FU", "AlcoholYN_FU",
               "HTSelfBP_FU", "DMSelfGluc_FU", "MetSyn_FU", "Lipidlowering_FU",
               "Antidepressants_FU", "Psychotropics_FU",
               "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU",
               "OralHealth_FU"),
    nose = c("Age_FU", "Sex", "BMI_FU", "Smoking_FU", "AlcoholYN_FU",
             "MetSyn_FU",
             "ToothBrushing_FU", "TongueBrushing_FU", "Mouthwash_FU")
)

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
## Returns a one-row summary tibble for the site-level pair summary.
run_da_pair <- function(ps_site, site_name, group1, group2, sig_confounders, outdir) {
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
        select(EthnicityTotal, all_of(sig_confounders))

    ## Fixed effects: ethnicity + significant confounders
    fixed_effects <- c("EthnicityTotal", sig_confounders)

    ## Output directory
    maaslin_outdir <- file.path(outdir, "maaslin2",
                                 paste0("maaslin2_16s_", site_name, "_", pair_tag))

    ## MaAsLin2's BH correction is computed independently within this model
    ## (across taxa). Unlike script 7's pairwise PERMANOVA - where every pair
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

    ## Save ethnicity-specific results
    write_csv(res_eth,
              file.path(outdir, "maaslin2",
                        paste0("maaslin2_ethnicity_results_16s_", site_name,
                               "_", pair_tag, ".csv")))

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

    ## =========================================================================
    ## Step 3: Visualization
    ## =========================================================================

    ## ---- Volcano plot ----
    res_eth <- res_eth |>
        mutate(
            sig = case_when(
                qval < 0.05 ~ "q < 0.05",
                qval < 0.25 ~ "q < 0.25",
                TRUE ~ "NS"
            ),
            ## Use the cleaned taxonomy label (built in 1b_datacleaning_biome.R)
            label = case_when(
                !is.na(Tax) & Tax != "" ~ Tax,
                TRUE ~ feature
            )
        ) |>
        ## Make labels unique by appending ASV ID for duplicates
        group_by(label) |>
        mutate(label = if (n() > 1) paste0(label, " (", feature, ")") else label) |>
        ungroup()

    volcano_colours <- c("q < 0.05" = "#E31A1C", "q < 0.25" = "#FF7F00",
                         "NS" = "grey60")

    p_volcano <- ggplot(res_eth, aes(x = coef, y = -log10(qval), colour = sig)) +
        geom_point(alpha = 0.7, size = 1.5) +
        geom_hline(yintercept = -log10(0.25), linetype = "dashed", colour = "grey40") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
        scale_colour_manual(values = volcano_colours) +
        labs(x = paste0("Coefficient (", group2, " vs ", group1, ")"),
             y = expression(-log[10](q-value)),
             colour = "Significance",
             title = paste0("Differential abundance - 16S ", site_name,
                            " (", g1_abbr, " vs ", g2_abbr, ")"),
             subtitle = paste0("MaAsLin2, adjusted for:\n",
                               format_covariates(sig_confounders))) +
        theme_Publication() +
        theme(legend.position = "right",
              plot.subtitle = element_text(size = rel(0.5)))

    ## Add labels for top significant taxa
    top_taxa <- res_eth |>
        filter(qval < 0.25) |>
        slice_min(qval, n = 15)

    if (nrow(top_taxa) > 0) {
        p_volcano <- p_volcano +
            ggrepel::geom_text_repel(
                data = top_taxa,
                aes(label = label),
                size = 2.5, max.overlaps = 20,
                show.legend = FALSE
            )
    }

    ggsave(file.path(outdir, "volcano",
                     paste0("volcano_16s_", site_name, "_", pair_tag, ".pdf")),
           plot = p_volcano, width = 8, height = 6)

    ## ---- Heatmap of top significant taxa ----
    if (n_sig > 0) {
        top_for_heatmap <- res_eth |>
            filter(qval < 0.25) |>
            slice_min(qval, n = min(30, n_sig))

        ## Get relative abundance for these taxa
        ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
        abund_mat <- as.data.frame(otu_table(ps_rel))
        if (taxa_are_rows(ps_rel)) {
            abund_mat <- abund_mat[top_for_heatmap$feature, , drop = FALSE]
        } else {
            abund_mat <- t(abund_mat)[top_for_heatmap$feature, , drop = FALSE]
        }

        ## Compute mean abundance per ethnicity. Both sapply and vapply
        ## silently simplify a length-1 per-group result down to a bare
        ## vector (losing the row dimension) when only one taxon is
        ## significant, so build the matrix explicitly via cbind instead,
        ## which always preserves it regardless of row count.
        mean_abund <- sapply(levels(meta$EthnicityTotal), function(eth) {
            samples <- rownames(meta[meta$EthnicityTotal == eth, ])
            rowMeans(abund_mat[, samples, drop = FALSE])
        }, simplify = FALSE) |> do.call(cbind, args = _)

        ## Use genus names for row labels
        rownames(mean_abund) <- top_for_heatmap$label

        ## Annotation: direction of effect, coloured by each group's
        ## canonical ethnicity colour so annotation colours stay consistent
        ## with the boxplots/PCoA across the whole analysis
        row_annotation <- data.frame(
            Direction = ifelse(top_for_heatmap$coef > 0,
                               paste("Enriched in", group2),
                               paste("Enriched in", group1)),
            row.names = top_for_heatmap$label
        )
        ann_colours <- list(Direction = setNames(
            c(eth_colours[[group2]], eth_colours[[group1]]),
            c(paste("Enriched in", group2), paste("Enriched in", group1))
        ))

        pdf(file.path(outdir, "heatmap",
                      paste0("heatmap_16s_", site_name, "_", pair_tag, ".pdf")),
            width = 9, height = max(4, nrow(mean_abund) * 0.3 + 2))
        pheatmap(log10(mean_abund + 1e-6),
                 cluster_cols = FALSE,
                 ## hclust needs >= 2 rows; a single significant taxon can't
                 ## be clustered against anything
                 cluster_rows = nrow(mean_abund) > 1,
                 annotation_row = row_annotation,
                 annotation_colors = ann_colours,
                 main = paste0("DA taxa - 16S ", site_name, " (", g1_abbr, " vs ",
                               g2_abbr, ")\n(q < 0.25, n = ", nrow(mean_abund), ")"),
                 fontsize_row = 8)
        dev.off()
    }

    ## ---- Forest plot: unadjusted vs adjusted, ordered by adjusted beta ----
    ## Taxa significant (q < 0.05) in the adjusted model; unadjusted estimates
    ## are shown alongside for comparison.
    sig_features <- res_eth |> filter(qval < 0.05) |> pull(feature)

    if (length(sig_features) > 0) {
        sig_taxa <- res_eth |> filter(feature %in% sig_features)

        forest_df <- bind_rows(
            sig_taxa |>
                select(feature, label, coef, stderr) |>
                mutate(model = "Adjusted"),
            res_eth_unadj |>
                filter(feature %in% sig_features) |>
                left_join(sig_taxa |> select(feature, label), by = "feature") |>
                select(feature, label, coef, stderr) |>
                mutate(model = "Unadjusted")
        ) |>
            mutate(
                conf.low  = coef - 1.96 * stderr,
                conf.high = coef + 1.96 * stderr
            )

        ## Order taxa by adjusted coefficient, most enriched in group2 first
        ## (highest coef), most depleted last (lowest coef).
        taxon_order <- sig_taxa |> arrange(desc(coef)) |> pull(label)
        forest_df <- forest_df |>
            mutate(label = factor(label, levels = taxon_order),
                   model = factor(model, levels = c("Unadjusted", "Adjusted")))

        model_colours <- c("Unadjusted" = "grey50", "Adjusted" = "#E31A1C")

        ## Paginate: a single page is capped at 50 inches by ggsave/pdf, so
        ## split into pages of at most 70 taxa (~23 in tall) when the
        ## significant set is large, keeping every taxon in the output.
        max_per_page <- 70
        n_total <- length(taxon_order)
        n_pages <- ceiling(n_total / max_per_page)

        pdf(file.path(outdir, "forestplot",
                      paste0("forestplot_16s_", site_name, "_", pair_tag, ".pdf")),
            width = 8, height = max(4, min(max_per_page, n_total) * 0.3 + 2))

        for (page in seq_len(n_pages)) {
            idx_start <- (page - 1) * max_per_page + 1
            idx_end   <- min(page * max_per_page, n_total)
            ## taxon_order is highest-to-lowest; reverse each page's slice so
            ## that after coord_flip() the highest value still lands at the
            ## top of the page (coord_flip puts the LAST factor level on top)
            page_taxa <- rev(taxon_order[idx_start:idx_end])

            page_df <- forest_df |>
                filter(label %in% page_taxa) |>
                mutate(label = factor(label, levels = page_taxa))

            p_forest <- ggplot(page_df, aes(x = label, y = coef, colour = model)) +
                geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
                geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                                position = position_dodge(width = 0.6),
                                size = 0.4, fatten = 1.5) +
                coord_flip() +
                scale_colour_manual(values = model_colours) +
                labs(x = NULL,
                     y = paste0("Coefficient (", group2, " vs ", group1, ", 95% CI)"),
                     colour = "Model",
                     title = paste0("Forest plot differential abundant ASVs ",
                                    site_name, " (", g1_abbr, " vs ", g2_abbr,
                                    ") (page ", page, "/", n_pages, ")"),
                     subtitle = paste0("q < 0.05 in adjusted model (n = ",
                                       n_total, " total)\nAdjusted for: ",
                                       format_covariates(sig_confounders))) +
                theme_Publication() +
                theme(legend.position = "right",
                      axis.text.y = element_text(size = rel(0.7)),
                      plot.subtitle = element_text(size = rel(0.55)))

            print(p_forest)
        }
        dev.off()
    }

    ## ---- Boxplots of top differentially abundant taxa ----
    if (n_strict > 0) {
        top_box <- res_eth |>
            filter(qval < 0.05) |>
            slice_min(qval, n = min(12, n_strict))

        ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
        abund_df <- as.data.frame(otu_table(ps_rel))
        if (taxa_are_rows(ps_rel)) {
            abund_df <- as.data.frame(t(abund_df))
        }

        ## Select top taxa and pivot
        box_df <- abund_df[, top_box$feature, drop = FALSE] |>
            rownames_to_column("sample_id") |>
            pivot_longer(-sample_id, names_to = "feature", values_to = "rel_abund") |>
            left_join(
                meta |> rownames_to_column("sample_id") |>
                    select(sample_id, EthnicityTotal),
                by = "sample_id"
            ) |>
            left_join(top_box |> select(feature, label, qval), by = "feature") |>
            mutate(label = paste0(label, "\n(q = ",
                                  formatC(qval, format = "e", digits = 1), ")"))

        ## Small pseudocount so zero-abundance samples remain visible on log scale
        ggplot(box_df, aes(x = EthnicityTotal, y = rel_abund + 1e-6, fill = EthnicityTotal)) +
            geom_boxplot(outlier.shape = 21, outlier.size = 0.5, alpha = 0.7) +
            facet_wrap(~ label, scales = "free_y") +
            scale_fill_manual(values = eth_colours) +
            scale_y_log10() +
            labs(x = NULL, y = "Relative abundance (log10 scale)", fill = "Ethnicity",
                 title = paste0("Top DA taxa - 16S ", site_name, " (", g1_abbr, " vs ",
                                g2_abbr, ", q < 0.05)")) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1),
                  strip.text = element_text(size = rel(0.6)))
        ggsave(file.path(outdir, "boxplots",
                        paste0("boxplots_top_da_16s_", site_name, "_", pair_tag, ".pdf")),
               width = 10, height = 8)
    }

    cat("Completed:", site_name, pair_tag, "\n")

    tibble(site = site_name, group1 = group1, group2 = group2,
           n_group1 = n_group1, n_group2 = n_group2,
           n_total = n_group1 + n_group2)
}

## ---- Analysis loop over sites ----
sites <- list(
    throat = readRDS("data/processed/ps_throat_rarefied.RDS"),
    nose   = readRDS("data/processed/ps_nose_rarefied.RDS")
)

## Nose has substantially fewer samples per group than throat (its largest
## non-Dutch group is Moroccan at 238 vs throat's 415+), so it uses its own,
## less stringent qualification threshold rather than throat's >400.
min_group_n <- c(throat = 400, nose = 150)
summary_rows <- list()

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
    write_csv(confounder_results,
              file.path(outdir, "confounder_assessment",
                        paste0("confounder_assessment_16s_", site_name, ".csv")))

    sig_confounders <- confounder_results |>
        filter(p.value < 0.05) |>
        pull(covariate)

    cat(site_name, "- Significant confounders (all groups):",
        paste(sig_confounders, collapse = ", "), "\n")

    for (pair in pairs) {
        summary_rows[[length(summary_rows) + 1]] <-
            run_da_pair(ps, site_name, pair[1], pair[2], sig_confounders, outdir)
    }

    cat("Completed:", site_name, "(", length(pairs), "pairs )\n\n")
}

write_csv(bind_rows(summary_rows), file.path(outdir, "summary_pairs_16s.csv"))
