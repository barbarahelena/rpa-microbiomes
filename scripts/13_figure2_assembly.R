## Figure 2 assembly: 16S throat/nose, stratified by ethnicity
## Panels:
##   A/B - Pairwise PERMANOVA (Weighted UniFrac, unadjusted) - nose/throat
##   C/D - Covariate effects on beta diversity (both distance metrics) - nose/throat
## Also writes a supplement (figure2_supp_braycurtis_pairwise.pdf/.png) with
## Bray-Curtis pairwise PERMANOVA heatmaps: A/B unadjusted (nose/throat),
## C/D adjusted for significant covariates (nose/throat).
## Rebuilt from the .rds cache written by 7a_beta_diversity_16s_compute.R, so
## this script never repeats a permutation test - run 7a first
## (BETA_DIV_TEST_N=<n> for a fast test cache).

## Libraries
library(here)
library(tidyverse)
library(cowplot)
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

## Word-wraps a title/subtitle to a fixed character width so it fits inside
## the plot instead of running off the page (same helper as 7b's report script)
wrap_for_plot <- function(x, width = 50) {
    if (is.null(x)) return(list(text = NULL, n_lines = 0))
    wrapped <- str_wrap(x, width = width)
    list(text = wrapped, n_lines = lengths(regmatches(wrapped, gregexpr("\n", wrapped))) + 1)
}

## Pairwise PERMANOVA R2 heatmap (BH-adjusted significance stars) - same as
## pairwise_permanova_heatmap() in 7b_beta_diversity_16s_report.R
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
    subtitle_wrapped <- wrap_for_plot(subtitle, width = 62)

    ggplot(pairwise_mat_df, aes(x = group1, y = group2, fill = R2)) +
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
}

## Setup
setwd(here::here())
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

## Beta diversity cache must match whatever 7a_beta_diversity_16s_compute.R
## was run with (BETA_DIV_TEST_N for a fast test cache, unset for the real one)
test_n <- suppressWarnings(as.integer(Sys.getenv("BETA_DIV_TEST_N", "")))
beta_outdir <- if (!is.na(test_n)) "results/beta_diversity_test" else "results/beta_diversity"
if (!is.na(test_n)) cat("TEST MODE: reading beta diversity cache from", beta_outdir, "\n")

## Colours for the two distance metrics in the covariate-effect summary plot
dist_colours <- c("Bray-Curtis" = "#2a78d6", "Weighted UniFrac" = "#eb6834")

## Human-readable labels for the covariate-effect barplot
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

sites <- c("nose", "throat")

heatmap_panels <- list()
adj_heatmap_panels <- list()
barplot_panels <- list()
n_terms_by_site <- list()
legend_source <- NULL
caches <- list()

for (site_name in sites) {
    cache_path <- file.path(beta_outdir, "cache", paste0("beta_diversity_16s_", site_name, ".rds"))
    if (!file.exists(cache_path)) {
        stop("No beta diversity cache for '", site_name, "' at ", cache_path,
             " - run scripts/7a_beta_diversity_16s_compute.R first.")
    }
    cache <- readRDS(cache_path)
    caches[[site_name]] <- cache
    meta      <- cache$meta
    distances <- cache$distances
    blocks    <- cache$blocks
    groups_order <- levels(droplevels(meta$EthnicityTotal))

    ## ---- Panel A/B: pairwise PERMANOVA heatmap, Weighted UniFrac, unadjusted ----
    wu_block <- blocks[["Weighted UniFrac"]]
    heatmap_panels[[site_name]] <- pairwise_permanova_heatmap(
        wu_block$permanova_pairwise, groups_order,
        title = paste0("Pairwise differences - ", site_name),
        fill_limits = c(0, max(wu_block$permanova_pairwise$R2))
    )

    ## ---- Panel E/F: pairwise PERMANOVA heatmap, Weighted UniFrac, adjusted
    ## for that site's significant covariates (see pairwise_permanova_adjusted()
    ## in 7a_beta_diversity_16s_compute.R). No subtitle listing the covariates
    ## here - same fill scale as A/B (not its own max) so the shrinkage after
    ## adjustment is visible directly by colour, title-only like A/B.
    adj_heatmap_panels[[site_name]] <- pairwise_permanova_heatmap(
        wu_block$permanova_pairwise_adjusted, groups_order,
        title = paste0("Pairwise differences (adjusted) - ", site_name),
        fill_limits = c(0, max(wu_block$permanova_pairwise$R2))
    )

    ## ---- Panel C/D: covariate effects on beta diversity, both distances ----
    covariate_screen_all <- list()
    for (dist_name in names(distances)) {
        block <- blocks[[dist_name]]
        covariate_screen_all[[dist_name]] <- bind_rows(
            tibble(covariate = "EthnicityTotal", Df = block$permanova_eth$Df[1],
                   R2 = block$permanova_eth$R2[1], F_stat = block$permanova_eth$F[1],
                   p.value = block$permanova_eth[["Pr(>F)"]][1]),
            block$covariate_screen
        ) |> mutate(distance = dist_name)
    }
    covariate_effects <- bind_rows(covariate_screen_all)

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

    n_terms_by_site[[site_name]] <- n_distinct(covariate_effects$covariate_label)

    p_bar <- ggplot(covariate_effects, aes(x = covariate_label, y = R2, fill = distance,
                                            colour = distance, alpha = significant)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.6) +
        coord_flip() +
        scale_fill_manual(values = dist_colours) +
        scale_colour_manual(values = dist_colours, guide = "none") +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
        labs(x = NULL, y = expression(R^2), fill = "Distance metric",
             title = paste0("Covariate effects - ", site_name),
             caption = "Filled bars: p < 0.05 (PERMANOVA)") +
        theme_Publication()

    if (is.null(legend_source)) legend_source <- p_bar

    barplot_panels[[site_name]] <- p_bar + theme(legend.position = "none")

    cat("Completed:", site_name, "\n")
}

shared_legend <- get_legend(legend_source + theme(legend.position = "bottom"))

## ---- Assemble ----
## Covariate barplots carry far more categories than the heatmaps have
## ethnicity groups, so the middle row gets proportionally more height.
grid <- plot_grid(
    heatmap_panels[["nose"]], heatmap_panels[["throat"]],
    barplot_panels[["nose"]], barplot_panels[["throat"]],
    adj_heatmap_panels[["nose"]], adj_heatmap_panels[["throat"]],
    ncol = 2, labels = c("A", "B", "C", "D", "E", "F"),
    rel_heights = c(0.8, 1.6, 0.8)
)

fig2 <- plot_grid(grid, shared_legend, ncol = 1, rel_heights = c(1, 0.04))

ggsave("results/figures/figure2.pdf", plot = fig2, width = 14, height = 21.5)
ggsave("results/figures/figure2.png", plot = fig2, width = 14, height = 21.5, dpi = 300)

cat("Saved results/figures/figure2.pdf and .png\n")

## ---- Supplement: Bray-Curtis pairwise PERMANOVA, unadjusted vs adjusted ----
## Reuses the cache objects already loaded above - no recompute. Unadjusted
## and adjusted heatmaps for the same site share one colour scale (as in
## 7b's per-site heatmaps) so the two are visually comparable.
bc_pairwise_panels <- list()

for (site_name in sites) {
    cache <- caches[[site_name]]
    meta  <- cache$meta
    block <- cache$blocks[["Bray-Curtis"]]
    groups_order <- levels(droplevels(meta$EthnicityTotal))

    shared_max <- max(block$permanova_pairwise$R2, block$permanova_pairwise_adjusted$R2)

    adjusted_subtitle <- if (length(block$sig_covariates) > 0) {
        paste0("Adjusted for: ", paste(covariate_labels[block$sig_covariates], collapse = ", "))
    } else {
        "No significant covariates - same as unadjusted"
    }

    bc_pairwise_panels[[paste0(site_name, "_unadj")]] <- pairwise_permanova_heatmap(
        block$permanova_pairwise, groups_order,
        title = paste0("Bray-Curtis pairwise (unadjusted) - ", site_name),
        fill_limits = c(0, shared_max)
    )
    bc_pairwise_panels[[paste0(site_name, "_adj")]] <- pairwise_permanova_heatmap(
        block$permanova_pairwise_adjusted, groups_order,
        title = paste0("Bray-Curtis pairwise (adjusted) - ", site_name),
        subtitle = adjusted_subtitle,
        fill_limits = c(0, shared_max)
    )
}

bc_grid <- plot_grid(
    bc_pairwise_panels[["nose_unadj"]],   bc_pairwise_panels[["throat_unadj"]],
    bc_pairwise_panels[["nose_adj"]],     bc_pairwise_panels[["throat_adj"]],
    ncol = 2, labels = c("A", "B", "C", "D")
)

## Fixed height (not dynamically sized to the adjusted panels' subtitle,
## unlike 7b's per-panel plots) - generous enough for the wrapped covariate
## list subtitle on the adjusted row
ggsave("results/figures/figure2_supp_braycurtis_pairwise.pdf", plot = bc_grid,
       width = 12, height = 13)
ggsave("results/figures/figure2_supp_braycurtis_pairwise.png", plot = bc_grid,
       width = 12, height = 13, dpi = 300)

cat("Saved results/figures/figure2_supp_braycurtis_pairwise.pdf and .png\n")

## ---- Supplement: ethnicity attenuation by covariate, nose vs throat ----
## Reuses the cache objects already loaded above - no recompute. Same
## summary plot as 7b_beta_diversity_16s_report.R's per-site attenuation
## plot (abs_reduction = how much ethnicity's PERMANOVA R2 drops once that
## one covariate is adjusted for), just nose and throat assembled side by
## side instead of as two separate per-site files.
attenuation_panels <- list()
attenuation_legend_source <- NULL

for (site_name in sites) {
    blocks <- caches[[site_name]]$blocks

    attenuation_effects <- map(names(blocks), function(dist_name) {
        blocks[[dist_name]]$ethnicity_attenuation |> mutate(distance = dist_name)
    }) |> bind_rows()

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

    p_atten <- ggplot(attenuation_effects, aes(x = covariate_label, y = abs_reduction,
                                                fill = distance, colour = distance,
                                                alpha = significant)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.6) +
        coord_flip() +
        scale_fill_manual(values = dist_colours) +
        scale_colour_manual(values = dist_colours, guide = "none") +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
        labs(x = NULL, y = expression("Reduction in ethnicity" ~ R^2),
             fill = "Distance metric",
             title = paste0("Ethnicity attenuation - ", site_name),
             caption = "Filled bars: p < 0.05 for ethnicity net of covariate (PERMANOVA)") +
        theme_Publication()

    if (is.null(attenuation_legend_source)) attenuation_legend_source <- p_atten

    attenuation_panels[[site_name]] <- p_atten + theme(legend.position = "none")
}

attenuation_legend <- get_legend(attenuation_legend_source + theme(legend.position = "bottom"))

n_atten_terms <- max(vapply(attenuation_panels,
                            function(p) n_distinct(p$data$covariate_label), integer(1)))

atten_grid <- plot_grid(
    attenuation_panels[["nose"]], attenuation_panels[["throat"]],
    ncol = 2, labels = c("A", "B")
)
atten_fig <- plot_grid(atten_grid, attenuation_legend, ncol = 1, rel_heights = c(1, 0.06))

ggsave("results/figures/figure2_supp_ethnicity_attenuation.pdf", plot = atten_fig,
       width = 14, height = max(6, 0.3 * n_atten_terms + 2))
ggsave("results/figures/figure2_supp_ethnicity_attenuation.png", plot = atten_fig,
       width = 14, height = max(6, 0.3 * n_atten_terms + 2), dpi = 300)

cat("Saved results/figures/figure2_supp_ethnicity_attenuation.pdf and .png\n")
