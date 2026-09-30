## Report ASV differential abundance from the matching computation cache.
## Export existing ethnicity-effect tables and volcano/heatmap/forest/boxplots.
## Includes the existing UpSet overlap report.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(pheatmap)

## Functions
source(here::here("scripts", "lib", "plot_style.R"))

## Setup
setwd(here::here())

## Read the matching full-data or DIFF_AB_TEST_N computation cache.
test_n <- suppressWarnings(as.integer(Sys.getenv("DIFF_AB_TEST_N", "")))
outdir <- if (!is.na(test_n)) "results/differential_abundance_test" else "results/differential_abundance"
if (!is.na(test_n)) cat("TEST MODE: reading saved results from", outdir, "\n")

for (sub in c("confounder_assessment", "maaslin2", "volcano", "heatmap", "forestplot", "boxplots")) {
    dir.create(file.path(outdir, sub), recursive = TRUE, showWarnings = FALSE)
}

## Define ethnicity colours
eth_colours <- ethnicity_colours

format_covariates <- function(covs, width = 70) {
    covs |>
        gsub("_(FU|BA)$", "", x = _) |>
        paste(collapse = ", ") |>
        str_wrap(width = width)
}

report_da_pair <- function(cache, outdir) {
    ps <- cache$ps
    meta <- cache$meta
    site_name <- cache$site_name
    group1 <- cache$group1
    group2 <- cache$group2
    sig_confounders <- cache$sig_confounders
    pair_tag <- cache$pair_tag
    g1_abbr <- cache$g1_abbr
    g2_abbr <- cache$g2_abbr
    n_group1 <- cache$n_group1
    n_group2 <- cache$n_group2
    res_eth <- cache$res_eth
    res_eth_unadj <- cache$res_eth_unadj
    n_sig <- cache$n_sig
    n_strict <- cache$n_strict
    volcano_seed <- cache$volcano_seed
    write_csv(res_eth,
              file.path(outdir, "maaslin2",
                        paste0("maaslin2_ethnicity_results_16s_", site_name,
                               "_", pair_tag, ".csv")))

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
            ## Use the cleaned taxonomy label (built in 02_clean_microbiome.R)
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
                seed = volcano_seed,
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
                     ## Just the comparison, in the same group2-vs-group1
                     ## order as the coefficient axis label above (titles
                     ## used to read group1-vs-group2, the reverse) - the
                     ## repeated "Forest plot differential abundant ASVs
                     ## <site>" boilerplate is dropped since these panels are
                     ## meant to sit side by side across all pairs at a site
                     title = if (n_pages > 1) {
                         paste0(g2_abbr, " vs ", g1_abbr,
                                " (page ", page, "/", n_pages, ")")
                     } else {
                         paste0(g2_abbr, " vs ", g1_abbr)
                     },
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
        ## Larger canvas than the panel count strictly needs, so the (often
        ## long) two-line taxon + q-value strip labels have room at the same
        ## font size instead of being truncated/overlapping.
        ggsave(file.path(outdir, "boxplots",
                        paste0("boxplots_top_da_16s_", site_name, "_", pair_tag, ".pdf")),
               width = 16, height = 11)
    }

}

cache_path <- file.path(outdir, "cache", "da_asv.rds")
if (!file.exists(cache_path)) stop("Missing cache: ", cache_path, "; run pixi run diffabund-16s-compute first.")
cache <- readRDS(cache_path)
for (site_name in names(cache$confounders)) {
    confounder_results <- cache$confounders[[site_name]]
    write_csv(confounder_results,
              file.path(outdir, "confounder_assessment",
                        paste0("confounder_assessment_16s_", site_name, ".csv")))

}
write_csv(cache$summary, file.path(outdir, "summary_pairs_16s.csv"))
for (pair_cache in cache$pairs) report_da_pair(pair_cache, outdir)

## Overlap of significant ASVs across pairs (no model fitting).
## Overlap of significant differential-abundance taxa across ethnicity pairs
## (16S throat and nose)
##
## Reads the per-pair MaAsLin2 results written above (maaslin2_ethnicity_results_16s_<site>_
## <pair>.csv) and builds an UpSet plot per site showing which ASVs are
## significant (q < 0.05) in which pairwise ethnicity comparisons, plus a CSV
## of the ASVs shared by 2+ pairs for follow-up. Doesn't refit anything, so
## it's cheap to re-run while iterating on the plot.

## Libraries
library(here)
library(tidyverse)
library(UpSetR)

## Setup
setwd(here::here())

## Use the same DIFF_AB_TEST_N output directory, so this can be pointed
## at a test-mode run's output for a quick check.
test_n <- suppressWarnings(as.integer(Sys.getenv("DIFF_AB_TEST_N", "")))
maaslin_dir <- if (!is.na(test_n)) "results/differential_abundance_test/maaslin2" else "results/differential_abundance/maaslin2"
outdir <- if (!is.na(test_n)) "results/differential_abundance_test/upset" else "results/differential_abundance/upset"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Significance threshold for set membership - matches the "strict" q < 0.05
## cutoff used for the forest plots and top-taxa boxplots above.
SIG_THRESHOLD <- 0.05

sites <- c("throat", "nose")

for (site_name in sites) {
    files <- list.files(
        maaslin_dir,
        pattern = paste0("^maaslin2_ethnicity_results_16s_", site_name, "_.+\\.csv$"),
        full.names = TRUE
    )

    if (length(files) == 0) {
        cat(site_name, "- no MaAsLin2 pair results found in", maaslin_dir,
            "- no pair results to report. Skipping.\n\n")
        next
    }

    ## pair_tag (e.g. "Dutch_vs_SAS") is everything after the site name
    pair_tags <- files |>
        basename() |>
        str_remove(paste0("^maaslin2_ethnicity_results_16s_", site_name, "_")) |>
        str_remove("\\.csv$")

    results <- map(files, read_csv, show_col_types = FALSE) |>
        set_names(pair_tags)

    ## Taxonomy lookup (feature -> readable label), same for every pair at a
    ## site since they share the same ASV table
    tax_lookup <- bind_rows(results) |>
        distinct(feature, Tax)

    ## Named list of significant ASVs per pair, for UpSetR::fromList()
    sig_sets <- map(results, ~ .x |> filter(qval < SIG_THRESHOLD) |> pull(feature))
    sig_sets <- sig_sets[lengths(sig_sets) > 0]

    if (length(sig_sets) < 2) {
        cat(site_name, "- fewer than 2 pairs have any significant taxa",
            "(q <", SIG_THRESHOLD, ") - skipping UpSet plot.\n\n")
        next
    }

    ## ---- UpSet plot ----
    membership <- UpSetR::fromList(sig_sets)
    pdf(file.path(outdir, paste0("upset_diffabund_16s_", site_name, ".pdf")),
        width = 9, height = 6, onefile = FALSE)
    print(UpSetR::upset(
        membership,
        sets = names(sig_sets),
        nintersects = 30,
        order.by = "freq",
        mainbar.y.label = paste0("Shared significant ASVs (q < ", SIG_THRESHOLD, ")"),
        sets.x.label = "Significant ASVs per pair",
        text.scale = 1.1
    ))
    dev.off()

    ## ---- Overlap detail: ASVs significant in 2+ pairs, with taxonomy ----
    ## fromList() doesn't keep the original element names as rownames - it
    ## just returns 1:N in the order of unique(unlist(sig_sets)), so rebuild
    ## that same vector to recover which feature each row corresponds to.
    membership_df <- membership |>
        mutate(feature = unique(unlist(sig_sets))) |>
        rowwise() |>
        mutate(n_pairs = sum(c_across(all_of(names(sig_sets)))),
               pairs = paste(names(sig_sets)[c_across(all_of(names(sig_sets))) == 1],
                            collapse = "; ")) |>
        ungroup() |>
        filter(n_pairs >= 2) |>
        left_join(tax_lookup, by = "feature") |>
        select(feature, Tax, n_pairs, pairs) |>
        arrange(desc(n_pairs))

    write_csv(membership_df,
              file.path(outdir, paste0("upset_overlap_16s_", site_name, ".csv")))

    cat("Completed:", site_name, "-", length(sig_sets), "pairs with signal,",
        nrow(membership_df), "ASVs shared by 2+ pairs\n\n")
}
