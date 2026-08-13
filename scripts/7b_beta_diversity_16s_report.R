## Beta diversity analysis: 16S microbiome (throat and nose) - REPORT STEP
## Reads the .rds cache written by 7a_beta_diversity_16s_compute.R and builds
## every plot and table (PCoA, betadisper, PERMANOVA tables/heatmaps,
## covariate screen, ethnicity attenuation) without repeating any permutation
## test. Run 5a first (or via the beta-16s pixi task, which chains both).

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

## Test mode must match the outdir 7a_beta_diversity_16s_compute.R wrote the
## cache to. Example: BETA_DIV_TEST_N=40 Rscript scripts/7b_beta_diversity_16s_report.R
test_n <- suppressWarnings(as.integer(Sys.getenv("BETA_DIV_TEST_N", "")))
outdir <- if (!is.na(test_n)) "results/beta_diversity_test" else "results/beta_diversity"
if (!is.na(test_n)) cat("TEST MODE: reading/writing", outdir, "\n")

for (sub in c("pcoa", "permanova", "covariate_screen", "betadisper")) {
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

## Colours for the two distance metrics in the covariate-effect summary plot
dist_colours <- c("Bray-Curtis" = "#2a78d6", "Weighted UniFrac" = "#eb6834")

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

## ---- Report loop over sites: throat and nose ----
for (site_name in c("throat", "nose")) {
    cache_path <- file.path(outdir, "cache", paste0("beta_diversity_16s_", site_name, ".rds"))
    if (!file.exists(cache_path)) {
        stop("No cache for '", site_name, "' at ", cache_path,
             " - run scripts/7a_beta_diversity_16s_compute.R first.")
    }
    cache <- readRDS(cache_path)
    meta      <- cache$meta
    n_samples <- cache$n_samples
    group_ns  <- cache$group_ns
    distances <- cache$distances
    pcoas     <- cache$pcoas
    blocks    <- cache$blocks

    cat("Groups (N>50) for", site_name, ":", paste(group_ns$EthnicityTotal, collapse = ", "), "\n")

    ## Collects ethnicity + covariate PERMANOVA R2/p from both distance
    ## metrics, for the covariate-effect summary plot built after this loop.
    covariate_screen_all <- list()

    ## Collects ethnicity-attenuation results from both distance metrics, for
    ## the attenuation summary plot built after this loop.
    ethnicity_attenuation_all <- list()

    for (dist_name in names(distances)) {
        dist_mat <- distances[[dist_name]]
        dist_label <- tolower(gsub("[- ]", "_", dist_name))
        pcoa <- pcoas[[dist_name]]

        block <- blocks[[dist_name]]
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
        ## in scripts/7a_beta_diversity_16s_compute.R)
        write_csv(ethnicity_attenuation,
                  paste0(outdir, "/covariate_screen/ethnicity_attenuation_", dist_label,
                         "_16s_", site_name, ".csv"))
        ethnicity_attenuation_all[[dist_name]] <- ethnicity_attenuation |>
            mutate(distance = dist_name)

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

    ## ---- Summary plot: ethnicity attenuation by covariate, across both
    ## distance metrics. abs_reduction is how much ethnicity's PERMANOVA R2
    ## drops once that one covariate is adjusted for - the covariate most
    ## responsible for confounding the ethnicity effect on beta diversity has
    ## the largest bar. ----
    attenuation_effects <- bind_rows(ethnicity_attenuation_all)

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
              paste0(outdir, "/permanova/ethnicity_attenuation_summary_16s_", site_name, ".csv"))

    n_atten_terms <- n_distinct(attenuation_effects$covariate_label)
    ggplot(attenuation_effects, aes(x = covariate_label, y = abs_reduction, fill = distance,
                                     colour = distance, alpha = significant)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.6) +
        coord_flip() +
        scale_fill_manual(values = dist_colours) +
        scale_colour_manual(values = dist_colours, guide = "none") +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
        labs(x = NULL, y = expression("Reduction in ethnicity" ~ R^2),
             fill = "Distance metric",
             title = paste0("Ethnicity attenuation by covariate - 16S ", site_name),
             caption = "Filled bars: p < 0.05 for ethnicity net of covariate (PERMANOVA)") +
        theme_Publication()
    ggsave(paste0(outdir, "/permanova/ethnicity_attenuation_summary_16s_", site_name, ".pdf"),
           width = 8, height = max(6, 0.3 * n_atten_terms + 2))

    cat("Finished site:", site_name, "-", n_samples, "samples\n\n")
}
