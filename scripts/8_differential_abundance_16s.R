## Differential abundance analysis: 16S microbiome (throat and nose)
## Ethnicity (Dutch vs South-Asian Surinamese)
## Step 1: Confounder assessment
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
dir.create("results/differential_abundance", recursive = TRUE, showWarnings = FALSE)

## Define ethnicity colours
eth_colours <- c("Dutch" = "#1F78B4", "South-Asian Surinamese" = "#E31A1C")

## Shorten covariate names for plot subtitles (drop trailing _FU/_BA suffix)
## and wrap long lists onto multiple lines instead of running off the page.
format_covariates <- function(covs, width = 70) {
    covs |>
        gsub("_(FU|BA)$", "", x = _) |>
        paste(collapse = ", ") |>
        str_wrap(width = width)
}

## Covariates to screen as potential confounders, per site.
## Site-specific lists follow the covariates that were significantly
## associated with beta diversity (Bray-Curtis and/or weighted UniFrac) for
## that site (see results/beta_diversity/covariate_screen_*_16s_<site>.csv).
## Antibiotics_FU is not screened here: antibiotic users are excluded from
## the analysis outright (see below) rather than adjusted for.
## Note: MigrationGen and ResidenceDuration_BA are excluded because they are
## structurally NA for all Dutch participants (migration-specific variables).
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

## ---- Analysis loop over sites ----
sites <- list(
    throat = readRDS("data/processed/ps_throat_rarefied.RDS"),
    nose   = readRDS("data/processed/ps_nose_rarefied.RDS")
)

for (site_name in names(sites)) {
    ps <- sites[[site_name]]
    covariates <- covariates_list[[site_name]]

    ## Filter to Dutch and South-Asian Surinamese only
    ps <- subset_samples(ps, EthnicityTotal %in% c("Dutch", "South-Asian Surinamese"))

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

    n_dutch <- sum(meta$EthnicityTotal == "Dutch")
    n_sas <- sum(meta$EthnicityTotal == "South-Asian Surinamese")

    ## =========================================================================
    ## Step 1: Confounder assessment
    ## Test association between each covariate and ethnicity
    ## =========================================================================
    confounder_results <- lapply(covariates, function(cov) {
        vals <- meta[[cov]]
        eth <- meta$EthnicityTotal
        cc <- !is.na(vals)
        vals <- vals[cc]
        eth <- eth[cc]

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
            res <- wilcox.test(vals ~ eth)
            return(tibble(covariate = cov, test = "wilcoxon",
                          statistic = res$statistic, p.value = res$p.value))
        }
    }) |> bind_rows()

    write_csv(confounder_results,
              paste0("results/differential_abundance/confounder_assessment_16s_",
                     site_name, ".csv"))

    ## Identify significant confounders (p < 0.05)
    sig_confounders <- confounder_results |>
        filter(p.value < 0.05) |>
        pull(covariate)

    cat("\n", site_name, "- Significant confounders:", paste(sig_confounders, collapse = ", "), "\n")

    ## =========================================================================
    ## Step 2: MaAsLin2 differential abundance
    ## =========================================================================

    ## Prepare input: taxa as columns, samples as rows
    counts_df <- as.data.frame(otu_table(ps))
    ## Ensure samples are rows
    if (taxa_are_rows(ps)) {
        counts_df <- as.data.frame(t(counts_df))
    }

    ## Prepare metadata for MaAsLin2
    meta_maaslin <- meta |>
        select(EthnicityTotal, all_of(sig_confounders))

    ## Fixed effects: ethnicity + significant confounders
    fixed_effects <- c("EthnicityTotal", sig_confounders)

    ## Output directory
    maaslin_outdir <- paste0("results/differential_abundance/maaslin2_16s_",
                             site_name)

    ## Run MaAsLin2
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
        reference      = c("EthnicityTotal,Dutch"),
        plot_heatmap   = FALSE,
        plot_scatter   = FALSE,
        max_significance = 0.25
    )

    ## Read results and filter to ethnicity effect
    res_all <- read_tsv(file.path(maaslin_outdir, "all_results.tsv"),
                        show_col_types = FALSE)
    res_eth <- res_all |>
        filter(metadata == "EthnicityTotal",
               value == "South-Asian Surinamese")

    ## Get taxonomy lookup
    tax <- as.data.frame(tax_table(ps)) |>
        rownames_to_column("feature")
    res_eth <- res_eth |>
        left_join(tax, by = "feature")

    ## Save ethnicity-specific results
    write_csv(res_eth,
              paste0("results/differential_abundance/maaslin2_ethnicity_results_16s_",
                     site_name, ".csv"))

    ## ---- Unadjusted model: ethnicity only, no confounders ----
    ## Fitted for comparison in the forest plot below.
    meta_maaslin_unadj <- meta |> select(EthnicityTotal)

    maaslin_outdir_unadj <- paste0("results/differential_abundance/maaslin2_16s_",
                                   site_name, "_unadjusted")

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
        reference      = c("EthnicityTotal,Dutch"),
        plot_heatmap   = FALSE,
        plot_scatter   = FALSE,
        max_significance = 0.25
    )

    res_eth_unadj <- read_tsv(file.path(maaslin_outdir_unadj, "all_results.tsv"),
                              show_col_types = FALSE) |>
        filter(metadata == "EthnicityTotal", value == "South-Asian Surinamese") |>
        select(feature, coef, stderr, pval, qval)

    n_sig <- sum(res_eth$qval < 0.25, na.rm = TRUE)
    n_strict <- sum(res_eth$qval < 0.05, na.rm = TRUE)
    cat(site_name, "- DA taxa (q < 0.25):", n_sig,
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
        labs(x = "Coefficient (South-Asian Surinamese vs Dutch)",
             y = expression(-log[10](q-value)),
             colour = "Significance",
             title = paste0("Differential abundance - 16S ", site_name),
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

    ggsave(paste0("results/differential_abundance/volcano_16s_", site_name, ".pdf"),
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
            abund_mat <- abund_mat[top_for_heatmap$feature, ]
        } else {
            abund_mat <- t(abund_mat)[top_for_heatmap$feature, ]
        }

        ## Compute mean abundance per ethnicity
        mean_abund <- sapply(levels(meta$EthnicityTotal), function(eth) {
            samples <- rownames(meta[meta$EthnicityTotal == eth, ])
            rowMeans(abund_mat[, samples, drop = FALSE])
        })

        ## Use genus names for row labels
        rownames(mean_abund) <- top_for_heatmap$label

        ## Annotation: direction of effect
        row_annotation <- data.frame(
            Direction = ifelse(top_for_heatmap$coef > 0, "Enriched in SAS", "Enriched in Dutch"),
            row.names = top_for_heatmap$label
        )
        ann_colours <- list(Direction = c("Enriched in SAS" = "#E31A1C",
                                          "Enriched in Dutch" = "#1F78B4"))

        pdf(paste0("results/differential_abundance/heatmap_16s_", site_name, ".pdf"),
            width = 6, height = max(4, nrow(mean_abund) * 0.3 + 2))
        pheatmap(log10(mean_abund + 1e-6),
                 cluster_cols = FALSE,
                 annotation_row = row_annotation,
                 annotation_colors = ann_colours,
                 main = paste0("DA taxa - 16S ", site_name,
                               " (q < 0.25, n = ", nrow(mean_abund), ")"),
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

        ## Order taxa by adjusted coefficient, most enriched in SAS first
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

        pdf(paste0("results/differential_abundance/forestplot_16s_", site_name, ".pdf"),
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
                     y = "Coefficient (South-Asian Surinamese vs Dutch, 95% CI)",
                     colour = "Model",
                     title = paste0("Forest plot differential abundant ASVs ", site_name,
                                    " (page ", page, "/", n_pages, ")"),
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
                 title = paste0("Top DA taxa - 16S ", site_name, " (q < 0.05)")) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 25, hjust = 1),
                  strip.text = element_text(size = rel(0.6)))
        ggsave(paste0("results/differential_abundance/boxplots_top_da_16s_",
                      site_name, ".pdf"),
               width = 10, height = 8)
    }

    cat("Completed:", site_name, "\n\n")
}
