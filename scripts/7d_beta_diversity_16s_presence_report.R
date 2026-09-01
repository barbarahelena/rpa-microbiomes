## Beta diversity analysis: 16S microbiome (throat and nose) - REPORT STEP
## PRESENCE/ABSENCE METRICS (unweighted UniFrac, binary Jaccard)
##
## Two jobs:
##
##   1. For the two presence/absence metrics, build the same per-metric plots
##      and tables 7b builds for the abundance-weighted ones (PCoA, betadisper,
##      PERMANOVA tables/heatmaps, covariate screen, ethnicity attenuation).
##      Reads the .rds cache written by
##      7c_beta_diversity_16s_presence_compute.R.
##
##   2. Build the four-metric comparison that motivates the whole analysis:
##      is ethnicity's effect on beta diversity carried by the abundant core
##      or by the rare tail? That needs the abundance-weighted numbers too, so
##      this script *also* reads the cache written by
##      7a_beta_diversity_16s_compute.R rather than recomputing them - the
##      weighted values here are then identical to the published ones in
##      results/beta_diversity/ by construction.
##
## No permutation test is repeated. Run 7a and 7c first (or via the
## beta-16s-presence pixi task, which chains what is needed).

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(grid)
library(ggthemes)

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

## Test mode must match the outdir the two compute scripts wrote their caches
## to - and note that BOTH must have been run in the same mode, since this
## script pairs their results.
## Example: BETA_DIV_TEST_N=40 Rscript scripts/7d_beta_diversity_16s_presence_report.R
test_n <- suppressWarnings(as.integer(Sys.getenv("BETA_DIV_TEST_N", "")))
outdir <- if (!is.na(test_n)) "results/beta_diversity_presence_test" else "results/beta_diversity_presence"
## Where 7a_beta_diversity_16s_compute.R put the abundance-weighted cache
weighted_outdir <- if (!is.na(test_n)) "results/beta_diversity_test" else "results/beta_diversity"
if (!is.na(test_n)) cat("TEST MODE: reading/writing", outdir, "\n")

for (sub in c("pcoa", "permanova", "covariate_screen", "betadisper", "comparison")) {
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

## Generic categorical palette for the supplementary covariate-coloured PCoA
## plots (same validated hues as eth_colours, applied by position since each
## covariate's factor levels differ)
cat_palette <- unname(eth_colours)

## Colours for the four distance metrics in the covariate-effect summary
## plots. The first two keep the exact hues 7b/13 use for them, so a reader
## comparing this figure against the published ones sees the same metric in
## the same colour. Hues are from the same validated palette as eth_colours.
dist_colours <- c(
    "Bray-Curtis"        = "#2a78d6",  # blue   - abundance,         taxonomic
    "Weighted UniFrac"   = "#eb6834",  # orange - abundance,         phylogenetic
    "Jaccard"            = "#33A02C",  # green  - presence/absence,  taxonomic
    "Unweighted UniFrac" = "#6A3D9A"   # purple - presence/absence,  phylogenetic
)

## Plot/legend order: the two abundance-weighted metrics first, then their
## presence/absence counterparts, so the comparison reads left-to-right.
metric_order <- c("Bray-Curtis", "Weighted UniFrac", "Jaccard", "Unweighted UniFrac")

## The metrics this script owns - the ones 7c computed. Everything else is
## read from 7a's cache purely for the comparison.
presence_metrics <- c("Unweighted UniFrac", "Jaccard")

## Which abundance-weighted metric each presence/absence metric is the
## counterpart of (same taxonomic-vs-phylogenetic family), for the
## weighted-vs-unweighted comparison table.
metric_counterpart <- c(
    "Jaccard"            = "Bray-Curtis",
    "Unweighted UniFrac" = "Weighted UniFrac"
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
    PM10_mean            = "PM10 (2013-2015 mean)",
    PM25_mean            = "PM2.5 (2013-2015 mean)",
    NO2_mean             = "NO2 (2014-2015 mean)",
    EC_mean              = "Soot/EC (2013-2015 mean)",
    Season               = "Collection season",
    EthnicityTotal       = "Ethnicity"
)

## Word-wraps a title/subtitle to a fixed character width so it fits inside
## the plot instead of running off the page, and reports how many lines the
## wrapped text ended up as (so callers can size the plot's height to it).
wrap_for_plot <- function(x, width = 50) {
    if (is.null(x)) return(list(text = NULL, n_lines = 0))
    wrapped <- str_wrap(x, width = width)
    list(text = wrapped, n_lines = lengths(regmatches(wrapped, gregexpr("\n", wrapped))) + 1)
}

## Pairwise PERMANOVA R2 heatmap (BH-adjusted significance stars). Each pair
## is drawn once - group1 always precedes group2 in groups_order (pairs come
## from combn() over an ordered group list) - so the tiles fill a single
## triangle instead of a mirrored square with every result shown twice.
## `fill_limits` lets an adjusted and unadjusted heatmap for the same site x
## metric share one colour scale, so a colour comparison between the two
## plots is fair rather than each auto-scaling to its own max R2.
pairwise_permanova_heatmap <- function(pairwise_df, groups_order, title, subtitle = NULL,
                                        fill_limits = c(0, NA)) {
    pairwise_mat_df <- pairwise_df |>
        select(group1, group2, R2, p.adj) |>
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

    title_wrapped <- wrap_for_plot(title, width = 45)
    ## Subtitle can carry a long comma-separated covariate list - render it
    ## smaller and wrap tighter than the title so lines stay inside a 6" plot
    subtitle_wrapped <- wrap_for_plot(subtitle, width = 62)

    p <- ggplot(pairwise_mat_df, aes(x = group1, y = group2, fill = R2)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = cell_label), size = 3, lineheight = 0.9) +
        scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", limits = fill_limits) +
        scale_x_discrete(drop = FALSE) +
        scale_y_discrete(drop = FALSE) +
        labs(x = NULL, y = NULL, fill = expression(R^2),
             title = title_wrapped$text, subtitle = subtitle_wrapped$text,
             caption = "BH-adjusted: * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right",
              panel.grid = element_blank(),
              plot.subtitle = element_text(size = rel(0.55), face = "italic"))

    ## Extra plot height (inches) the caller should add so the wrapped
    ## title/subtitle have room to breathe instead of crowding the tiles
    attr(p, "extra_height") <- 0.28 * (title_wrapped$n_lines - 1) +
        0.17 * subtitle_wrapped$n_lines
    p
}

## Loads one compute script's cache for one site, or stops with a message
## naming the script that produces it.
read_beta_cache <- function(dir, site_name, produced_by) {
    cache_path <- file.path(dir, "cache", paste0("beta_diversity_16s_", site_name, ".rds"))
    if (!file.exists(cache_path)) {
        stop("No cache for '", site_name, "' at ", cache_path,
             " - run scripts/", produced_by, " first",
             if (!is.na(test_n)) paste0(" (with BETA_DIV_TEST_N=", test_n, ")") else "",
             ".")
    }
    readRDS(cache_path)
}

## Collects the omnibus ethnicity result per site x metric, for the
## weighted-vs-unweighted markdown summary written after the site loop.
omnibus_all <- list()

## ---- Report loop over sites: throat and nose ----
for (site_name in c("throat", "nose")) {
    presence_cache <- read_beta_cache(outdir, site_name,
                                      "7c_beta_diversity_16s_presence_compute.R")
    weighted_cache <- read_beta_cache(weighted_outdir, site_name,
                                      "7a_beta_diversity_16s_compute.R")

    meta      <- presence_cache$meta
    n_samples <- presence_cache$n_samples
    group_ns  <- presence_cache$group_ns
    pcoas     <- presence_cache$pcoas
    blocks    <- presence_cache$blocks

    ## The comparison is only meaningful if both caches cover the same
    ## samples. They will whenever 7a and 7c were run in the same mode (the
    ## filtering code in 7c is identical to 7a's), so a mismatch here means
    ## the caches are from different runs - most likely one real and one
    ## BETA_DIV_TEST_N - and silently plotting them side by side would
    ## produce a wrong comparison.
    if (!identical(rownames(weighted_cache$meta), rownames(meta))) {
        stop("Cache mismatch for '", site_name, "': ", weighted_outdir, " has ",
             nrow(weighted_cache$meta), " samples and ", outdir, " has ", nrow(meta),
             " (or the same count in a different order). Re-run both ",
             "scripts/7a_beta_diversity_16s_compute.R and ",
             "scripts/7c_beta_diversity_16s_presence_compute.R in the same mode.")
    }

    ## All four metrics in one list for the cross-metric comparison plots.
    ## Per-metric outputs below are written only for the presence/absence
    ## ones - the weighted metrics already have theirs in results/beta_diversity/.
    all_blocks <- c(weighted_cache$blocks, blocks)[metric_order]

    cat("Groups (N>50) for", site_name, ":", paste(group_ns$EthnicityTotal, collapse = ", "), "\n")

    ## Collects ethnicity + covariate PERMANOVA R2/p from all four distance
    ## metrics, for the covariate-effect summary plot built after this loop.
    covariate_screen_all <- list()

    ## Collects ethnicity-attenuation results from all four distance metrics,
    ## for the attenuation summary plot built after this loop.
    ethnicity_attenuation_all <- list()

    for (dist_name in metric_order) {
        block <- all_blocks[[dist_name]]

        ## Stash effect sizes for the four-metric summary plots. Done for
        ## every metric, including the two read out of 7a's cache.
        covariate_screen_all[[dist_name]] <- bind_rows(
            tibble(covariate = "EthnicityTotal", Df = block$permanova_eth$Df[1],
                   R2 = block$permanova_eth$R2[1], F_stat = block$permanova_eth$F[1],
                   p.value = block$permanova_eth[["Pr(>F)"]][1]),
            block$covariate_screen
        ) |> mutate(distance = dist_name)

        ethnicity_attenuation_all[[dist_name]] <- block$ethnicity_attenuation |>
            mutate(distance = dist_name)

        ## Omnibus ethnicity R2/p, unadjusted and adjusted, for the markdown
        ## summary. permanova_full falls back to permanova_eth when a metric
        ## had no significant covariates, in which case the two agree.
        omnibus_all[[paste(site_name, dist_name)]] <- tibble(
            site            = site_name,
            distance        = dist_name,
            R2_unadjusted   = block$permanova_eth$R2[1],
            p_unadjusted    = block$permanova_eth[["Pr(>F)"]][1],
            R2_adjusted     = block$permanova_full["EthnicityTotal", "R2"],
            p_adjusted      = block$permanova_full["EthnicityTotal", "Pr(>F)"],
            n_sig_covariates = length(block$sig_covariates),
            betadisper_p    = block$betadisp_test$tab[["Pr(>F)"]][1]
        )

        ## Per-metric plots and tables are this script's own output only for
        ## the presence/absence metrics
        if (!dist_name %in% presence_metrics) next

        dist_label <- tolower(gsub("[- ]", "_", dist_name))
        pcoa <- pcoas[[dist_name]]

        permanova_eth               <- block$permanova_eth
        permanova_pairwise          <- block$permanova_pairwise
        permanova_pairwise_adjusted <- block$permanova_pairwise_adjusted
        covariate_screen            <- block$covariate_screen
        sig_covariates              <- block$sig_covariates
        ethnicity_attenuation       <- block$ethnicity_attenuation
        permanova_full              <- block$permanova_full
        betadisp                    <- block$betadisp
        betadisp_test               <- block$betadisp_test
        betadisp_pairwise           <- block$betadisp_pairwise

        permanova_label <- paste0(
            "PERMANOVA: R² = ", round(permanova_eth$R2[1], 3),
            ", p = ", format.pval(permanova_eth[["Pr(>F)"]][1], digits = 2, eps = 0.001)
        )

        ## ---- PCoA ordination ----
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

        ## PERMANOVA pairwise post-hoc, adjusted for significant covariates
        ## (same covariates as permanova_full's adjustment set for this site
        ## x distance combination; empty set means this equals the
        ## unadjusted table above)
        write_csv(permanova_pairwise_adjusted,
                  paste0(outdir, "/permanova/permanova_pairwise_adjusted_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## ---- Pairwise PERMANOVA heatmaps (R2, BH-adjusted significance) ----
        groups_order <- levels(droplevels(meta$EthnicityTotal))

        ## Shared colour-scale ceiling across the unadjusted and adjusted
        ## heatmaps, so a colour comparison between the two is fair instead
        ## of each plot auto-scaling to its own max R2.
        shared_fill_max <- max(permanova_pairwise$R2, permanova_pairwise_adjusted$R2)

        p_permanova_heat <- pairwise_permanova_heatmap(
            permanova_pairwise, groups_order,
            title = paste0("Pairwise PERMANOVA - ", dist_name, " - 16S ", site_name),
            fill_limits = c(0, shared_fill_max)
        )
        ggsave(paste0(outdir, "/permanova/permanova_pairwise_heatmap_", dist_label,
                      "_16s_", site_name, ".pdf"),
               plot = p_permanova_heat,
               width = 6, height = 5.5 + attr(p_permanova_heat, "extra_height"))

        adjusted_subtitle <- if (length(sig_covariates) > 0) {
            paste0("Adjusted for: ", paste(covariate_labels[sig_covariates], collapse = ", "))
        } else {
            "No significant covariates - same as unadjusted"
        }
        p_permanova_heat_adjusted <- pairwise_permanova_heatmap(
            permanova_pairwise_adjusted, groups_order,
            title = paste0("Pairwise PERMANOVA (adjusted) - ", dist_name, " - 16S ", site_name),
            subtitle = adjusted_subtitle,
            fill_limits = c(0, shared_fill_max)
        )
        ggsave(paste0(outdir, "/permanova/permanova_pairwise_heatmap_adjusted_", dist_label,
                      "_16s_", site_name, ".pdf"),
               plot = p_permanova_heat_adjusted,
               width = 6, height = 5.5 + attr(p_permanova_heat_adjusted, "extra_height"))

        ## Covariate screening
        write_csv(covariate_screen,
                  paste0(outdir, "/covariate_screen/covariate_screen_", dist_label,
                         "_16s_", site_name, ".csv"))

        ## Ethnicity attenuation: which covariates explain the most of
        ## ethnicity's effect on beta diversity (see ethnicity_attenuation()
        ## in scripts/7c_beta_diversity_16s_presence_compute.R)
        write_csv(ethnicity_attenuation,
                  paste0(outdir, "/covariate_screen/ethnicity_attenuation_", dist_label,
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

    ## ---- Comparison plot: ethnicity + covariate R2 across all four distance
    ## metrics. This is the plot the analysis exists for: ethnicity's bar being
    ## taller under Jaccard/unweighted UniFrac than under their abundance-
    ## weighted counterparts puts the effect in the rare tail, shorter puts it
    ## in the abundant core. ----
    covariate_effects <- bind_rows(covariate_screen_all) |>
        mutate(distance = factor(distance, levels = metric_order))

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
              paste0(outdir, "/comparison/permanova_summary_allmetrics_16s_", site_name, ".csv"))

    ## Significant bars are filled solid; non-significant bars show only
    ## their coloured outline (fill alpha = 0) - no text label needed, and
    ## it sidesteps ever having to align a star with a dodged bar.
    ## Four dodged bars per covariate rather than 7b's two, so the panel gets
    ## a little more width and height to keep the bars legible.
    n_terms <- n_distinct(covariate_effects$covariate_label)
    ggplot(covariate_effects, aes(x = covariate_label, y = R2, fill = distance,
                                   colour = distance, alpha = significant)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, linewidth = 0.4) +
        coord_flip() +
        scale_fill_manual(values = dist_colours) +
        scale_colour_manual(values = dist_colours, guide = "none") +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
        labs(x = NULL, y = expression(R^2), fill = "Distance metric",
             title = paste0("Covariate effects on beta diversity - 16S ", site_name),
             subtitle = "Abundance-weighted (blue, orange) vs presence/absence (green, purple) metrics",
             caption = "Filled bars: p < 0.05 (PERMANOVA)") +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_Publication() +
        theme(plot.subtitle = element_text(size = rel(0.7), hjust = 0.5))
    ggsave(paste0(outdir, "/comparison/permanova_summary_allmetrics_16s_", site_name, ".pdf"),
           width = 9, height = max(6, 0.4 * n_terms + 2))

    ## ---- Comparison plot: ethnicity attenuation by covariate, across all
    ## four distance metrics. abs_reduction is how much ethnicity's PERMANOVA
    ## R2 drops once that one covariate is adjusted for - the covariate most
    ## responsible for confounding the ethnicity effect on beta diversity has
    ## the largest bar. A covariate that confounds only the presence/absence
    ## metrics is acting on the rare tail specifically. ----
    attenuation_effects <- bind_rows(ethnicity_attenuation_all) |>
        mutate(distance = factor(distance, levels = metric_order))

    attenuation_order <- attenuation_effects |>
        group_by(covariate) |>
        summarise(max_reduction = max(abs_reduction), .groups = "drop") |>
        arrange(max_reduction) |>
        pull(covariate)
    attenuation_effects <- attenuation_effects |>
        mutate(covariate_label = covariate_labels[covariate],
               covariate_label = factor(covariate_label,
                                        levels = covariate_labels[attenuation_order]),
               significant = p_value < 0.05)

    write_csv(attenuation_effects,
              paste0(outdir, "/comparison/ethnicity_attenuation_allmetrics_16s_", site_name, ".csv"))

    n_atten_terms <- n_distinct(attenuation_effects$covariate_label)
    ggplot(attenuation_effects, aes(x = covariate_label, y = abs_reduction, fill = distance,
                                     colour = distance, alpha = significant)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, linewidth = 0.4) +
        coord_flip() +
        scale_fill_manual(values = dist_colours) +
        scale_colour_manual(values = dist_colours, guide = "none") +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
        labs(x = NULL, y = expression("Reduction in ethnicity" ~ R^2),
             fill = "Distance metric",
             title = paste0("Ethnicity attenuation by covariate - 16S ", site_name),
             subtitle = "Abundance-weighted (blue, orange) vs presence/absence (green, purple) metrics",
             caption = "Filled bars: p < 0.05 for ethnicity net of covariate (PERMANOVA)") +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_Publication() +
        theme(plot.subtitle = element_text(size = rel(0.7), hjust = 0.5))
    ggsave(paste0(outdir, "/comparison/ethnicity_attenuation_allmetrics_16s_", site_name, ".pdf"),
           width = 9, height = max(6, 0.4 * n_atten_terms + 2))

    cat("Finished site:", site_name, "-", n_samples, "samples\n\n")
}

## ---- Written comparison: weighted vs unweighted ethnicity effect ----
## Generated from the numbers just loaded rather than hand-written, so it
## cannot drift out of step with the data the way a manually maintained
## results summary does.
omnibus <- bind_rows(omnibus_all) |>
    mutate(distance = factor(distance, levels = metric_order)) |>
    arrange(site, distance)

fmt_pct <- function(x) sprintf("%.2f%%", 100 * x)
fmt_p   <- function(x) format.pval(x, digits = 2, eps = 0.001)

## Each presence/absence metric against its abundance-weighted counterpart
## (Jaccard vs Bray-Curtis, unweighted vs weighted UniFrac), so the ratio
## compares like with like - taxonomic against taxonomic, phylogenetic
## against phylogenetic.
comparison <- omnibus |>
    filter(distance %in% presence_metrics) |>
    transmute(
        site,
        family = if_else(distance == "Jaccard", "Taxonomic", "Phylogenetic"),
        presence_metric = as.character(distance),
        weighted_metric = unname(metric_counterpart[as.character(distance)])
    ) |>
    left_join(omnibus |> transmute(site, presence_metric = as.character(distance),
                                   R2_pa = R2_unadjusted, p_pa = p_unadjusted,
                                   R2_pa_adj = R2_adjusted, disp_p_pa = betadisper_p),
              by = c("site", "presence_metric")) |>
    left_join(omnibus |> transmute(site, weighted_metric = as.character(distance),
                                   R2_w = R2_unadjusted, p_w = p_unadjusted,
                                   R2_w_adj = R2_adjusted, disp_p_w = betadisper_p),
              by = c("site", "weighted_metric")) |>
    mutate(ratio = R2_pa / R2_w,
           ratio_adj = R2_pa_adj / R2_w_adj)

write_csv(comparison,
          file.path(outdir, "comparison", "weighted_vs_unweighted_ethnicity_R2.csv"))

## One sentence per site x family, stating which side of the abundance
## distribution the ethnicity effect sits on. The 1.1x / 0.9x band is a
## deliberately loose "no material difference" zone - PERMANOVA R2 from 999
## permutations on this many samples is not precise enough for a 5% gap
## between metrics to mean anything mechanistic.
interpret <- function(row) {
    verdict <- if (row$ratio >= 1.1) {
        paste0("**higher** under ", row$presence_metric, " than under ",
               row$weighted_metric, " (", sprintf("%.2fx", row$ratio),
               "), so the ethnicity signal is concentrated in the rare/low-abundance ",
               "tail: groups differ in *which* uncommon taxa they carry more than in ",
               "the relative abundances of the dominant ones")
    } else if (row$ratio <= 0.9) {
        paste0("**lower** under ", row$presence_metric, " than under ",
               row$weighted_metric, " (", sprintf("%.2fx", row$ratio),
               "), so the ethnicity signal is carried by shifts among the abundant ",
               "taxa rather than by turnover in the rare tail")
    } else {
        paste0("**comparable** under ", row$presence_metric, " and ",
               row$weighted_metric, " (", sprintf("%.2fx", row$ratio),
               "), so the effect is spread across the abundance distribution rather ",
               "than concentrated in either the core or the rare tail")
    }
    paste0("- **", row$site, ", ", tolower(row$family), " (",
           row$presence_metric, " vs ", row$weighted_metric, ")**: ethnicity R² is ",
           verdict, ". Unadjusted R² ", fmt_pct(row$R2_pa), " (p = ", fmt_p(row$p_pa),
           ") vs ", fmt_pct(row$R2_w), " (p = ", fmt_p(row$p_w),
           "); adjusted for that metric's significant covariates, ",
           fmt_pct(row$R2_pa_adj), " vs ", fmt_pct(row$R2_w_adj),
           " (", sprintf("%.2fx", row$ratio_adj), ").")
}

## Betadisper matters here on its own: dispersion that differs between groups
## only under a presence/absence metric means the groups differ in how
## variable their rare fraction is, which PERMANOVA alone cannot separate from
## a location shift.
interpret_dispersion <- function(row) {
    pa_sig <- row$disp_p_pa < 0.05
    w_sig  <- row$disp_p_w < 0.05
    note <- if (pa_sig && !w_sig) {
        paste0("dispersion differs between groups under ", row$presence_metric,
               " but not under ", row$weighted_metric,
               " - the groups differ in how variable their rare fraction is, so part of ",
               "the PERMANOVA signal here may be a spread difference rather than a ",
               "location shift")
    } else if (!pa_sig && w_sig) {
        paste0("dispersion differs under ", row$weighted_metric, " but not under ",
               row$presence_metric, " - the heterogeneity is in abundance, not in taxon membership")
    } else if (pa_sig && w_sig) {
        "dispersion differs between groups under both metrics, so the PERMANOVA result should not be read as a pure location shift"
    } else {
        "dispersion does not differ between groups under either metric, so the PERMANOVA result reads as a location shift"
    }
    paste0("- **", row$site, ", ", tolower(row$family), "**: ", note,
           " (betadisper permutest p = ", fmt_p(row$disp_p_pa), " vs ",
           fmt_p(row$disp_p_w), ").")
}

md <- c(
    "# Ethnicity effect on beta diversity: abundance-weighted vs presence/absence metrics",
    "",
    paste0("Generated by `scripts/7d_beta_diversity_16s_presence_report.R` on ",
           format(Sys.Date()), ". Do not edit by hand - re-run the script."),
    "",
    if (!is.na(test_n)) c(
        paste0("> **TEST MODE (BETA_DIV_TEST_N=", test_n, ").** Each ethnicity group was ",
               "capped at ", test_n, " samples. These numbers exist to check the pipeline ",
               "and are not interpretable."),
        ""
    ),
    "## Why",
    "",
    "Bray-Curtis and weighted UniFrac are both abundance-weighted, so they are dominated by",
    "the handful of abundant taxa and are the metrics least sensitive to a difference",
    "concentrated in the rare tail. Alpha diversity points at exactly such a difference:",
    "migrant groups carry more taxa at lower evenness than the Dutch group at throat.",
    "Jaccard and unweighted UniFrac weight every observed taxon equally, so comparing",
    "ethnicity's R² between the two families locates the effect on the abundance",
    "distribution. Each presence/absence metric is compared against its counterpart in the",
    "same family - Jaccard against Bray-Curtis (taxonomic), unweighted against weighted",
    "UniFrac (phylogenetic).",
    "",
    "## Omnibus ethnicity effect (PERMANOVA, 999 permutations)",
    "",
    "| Site | Metric | Weighting | Unadjusted R² | p | Adjusted R² | p | Sig. covariates |",
    "|---|---|---|---|---|---|---|---|",
    omnibus |>
        mutate(weighting = if_else(as.character(distance) %in% presence_metrics,
                                   "presence/absence", "abundance"),
               line = paste0("| ", site, " | ", distance, " | ", weighting, " | ",
                             fmt_pct(R2_unadjusted), " | ", fmt_p(p_unadjusted), " | ",
                             fmt_pct(R2_adjusted), " | ", fmt_p(p_adjusted), " | ",
                             n_sig_covariates, " |")) |>
        pull(line),
    "",
    "Adjusted R² is ethnicity net of that metric's own significant covariates, so the",
    "adjustment set differs between metrics and the adjusted columns are not a like-for-like",
    "comparison the way the unadjusted ones are.",
    "",
    "## Where does the effect sit?",
    "",
    vapply(seq_len(nrow(comparison)), function(i) interpret(comparison[i, ]), character(1)),
    "",
    "## Dispersion (betadisper)",
    "",
    vapply(seq_len(nrow(comparison)), function(i) interpret_dispersion(comparison[i, ]), character(1)),
    "",
    "## Outputs",
    "",
    paste0("- Per-metric tables and plots for the presence/absence metrics: `", outdir,
           "/{permanova,pcoa,betadisper,covariate_screen}/`"),
    paste0("- Four-metric covariate-effect and attenuation plots: `", outdir, "/comparison/`"),
    paste0("- The abundance-weighted numbers above are read unchanged from `",
           weighted_outdir, "/cache/`, so they match the published results in that folder."),
    ""
)

writeLines(md, file.path(outdir, "comparison", "weighted_vs_unweighted_summary.md"))
cat("Wrote", file.path(outdir, "comparison", "weighted_vs_unweighted_summary.md"), "\n")
