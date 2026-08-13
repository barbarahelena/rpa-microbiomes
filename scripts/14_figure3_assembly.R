## Figure 3 assembly: 16S throat/nose differential abundance (ethnicity)
## Main figure panels:
##   A/B - Pairwise DA taxa count heatmap (q < 0.05), nose/throat - summarizes
##         all 10/15 pairwise MaAsLin2 models at a glance
##   C     - Volcano plot, nose: Dutch vs African Surinamese (best nose signal)
##   D/E/F - Volcano plots, throat: Dutch vs South-Asian Surinamese (strongest
##           signal), Dutch vs Ghanaian (second-strongest), Ghanaian vs
##           Moroccan (non-Dutch pair with strong signal, shows migrant-
##           migrant divergence is real too)
## Also writes two supplementary grids (not part of the lettered main
## figure) with every pairwise volcano plot for a site, so the pairs not
## featured above are still available for the supplement.
## Volcano plots have no adjustment-set subtitle - state it in the figure
## legend instead (same adjustment set for every pair at a site, see
## results/differential_abundance/confounder_assessment/, printed below).
## Rebuilt from the per-pair MaAsLin2 CSVs already written by
## 9_differential_abundance_16s.R, so this script never refits anything -
## run 9 first (DIFF_AB_TEST_N=<n> for a fast test run, matched here).

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

## Shorten covariate names for plot subtitles (drop trailing _FU/_BA suffix)
## and wrap long lists onto multiple lines instead of running off the page.
## (same helper as 9_differential_abundance_16s.R)
format_covariates <- function(covs, width = 70) {
    covs |>
        gsub("_(FU|BA)$", "", x = _) |>
        paste(collapse = ", ") |>
        str_wrap(width = width)
}

## Short abbreviations for filenames (group names contain spaces/hyphens) -
## same lookup as 9_differential_abundance_16s.R, needed here to rebuild the
## pair_tag used in its output filenames.
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

eth_colours <- c(
    "Dutch"                  = "#1F78B4",
    "South-Asian Surinamese" = "#E31A1C",
    "African Surinamese"     = "#33A02C",
    "Javanese Surinamese"    = "#6A3D9A",
    "Other"                  = "#B15928",
    "Ghanaian"               = "#FF7F00",
    "Turkish"                = "#E7298A",
    "Moroccan"               = "#D4AC0D"
)

## Pairwise DA taxa count heatmap - same visual convention as the pairwise
## PERMANOVA heatmap in 13_figure2_assembly.R (upper-triangle tile grid,
## sequential fill), but cells hold a DA taxa count instead of R2.
pairwise_sig_heatmap <- function(pair_counts, groups_order, title, fill_limits = c(0, NA)) {
    pair_counts <- pair_counts |>
        mutate(group1 = factor(group1, levels = groups_order),
               group2 = factor(group2, levels = groups_order))

    ggplot(pair_counts, aes(x = group1, y = group2, fill = n_sig)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = n_sig), size = 3.5) +
        scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", limits = fill_limits) +
        scale_x_discrete(drop = FALSE) +
        scale_y_discrete(drop = FALSE) +
        labs(x = NULL, y = NULL, fill = "Differentially\nabundant microbes\n(q < 0.05)", title = title) +
        guides(fill = guide_colorbar(barheight = unit(1.8, "cm"), barwidth = unit(0.3, "cm"))) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right",
              panel.grid = element_blank())
}

## Volcano plot for a single pair, rebuilt from the MaAsLin2 CSV already
## written by 9_differential_abundance_16s.R (no refitting). Mirrors that
## script's volcano code, minus the file save and the adjustment-set
## subtitle (stated once in the figure legend instead - see
## site_confounder_label() below). label_n/title_site are lower/simpler
## for the all-pairs supplementary grids than for the main figure's
## highlighted panels.
build_volcano <- function(diffab_dir, site_name, group1, group2,
                           label_n = 10, title_site = TRUE) {
    pair_tag <- pair_name(group1, group2)
    g1_abbr <- eth_abbrev[[group1]]
    g2_abbr <- eth_abbrev[[group2]]

    res_eth <- read_csv(
        file.path(diffab_dir, "maaslin2",
                  paste0("maaslin2_ethnicity_results_16s_", site_name, "_", pair_tag, ".csv")),
        show_col_types = FALSE
    )

    res_eth <- res_eth |>
        mutate(
            sig = case_when(
                qval < 0.05 ~ "q < 0.05",
                TRUE ~ "NS"
            ),
            label = case_when(
                !is.na(Tax) & Tax != "" ~ Tax,
                TRUE ~ feature
            )
        ) |>
        group_by(label) |>
        mutate(label = if (n() > 1) paste0(label, " (", feature, ")") else label) |>
        ungroup()

    volcano_colours <- c("q < 0.05" = "#E31A1C", "NS" = "grey60")

    p <- ggplot(res_eth, aes(x = coef, y = -log10(qval), colour = sig)) +
        geom_point(alpha = 0.7, size = 1.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
        scale_colour_manual(values = volcano_colours) +
        labs(x = paste0("Coefficient (", group2, " vs ", group1, ")"),
             y = expression(-log[10](q-value)),
             colour = "Significance",
             title = if (title_site) paste0(site_name, ": ", g1_abbr, " vs ", g2_abbr)
                     else paste0(g1_abbr, " vs ", g2_abbr)) +
        theme_Publication() +
        theme(legend.position = "right")

    top_taxa <- res_eth |>
        filter(qval < 0.05) |>
        slice_min(qval, n = label_n)

    if (nrow(top_taxa) > 0) {
        p <- p +
            ggrepel::geom_text_repel(
                data = top_taxa,
                aes(label = label),
                size = 2.3, max.overlaps = 15,
                show.legend = FALSE
            )
    }

    p
}

## Significant-confounder set for a site, formatted for the figure legend
## text (same set adjusted for in every pairwise model at that site - see
## 9_differential_abundance_16s.R's assess_confounders()).
site_confounder_label <- function(diffab_dir, site_name) {
    sig_confounders <- read_csv(
        file.path(diffab_dir, "confounder_assessment",
                  paste0("confounder_assessment_16s_", site_name, ".csv")),
        show_col_types = FALSE
    ) |>
        filter(p.value < 0.05) |>
        pull(covariate)
    format_covariates(sig_confounders, width = 200)
}

## Grid of every pairwise volcano plot at a site (10 at nose, 15 at
## throat) - the supplementary counterpart to the main figure's 4
## highlighted pairs.
assemble_site_volcanoes <- function(diffab_dir, site_name, summary_pairs, ncol = 3, label_n = 6) {
    pairs <- summary_pairs |> filter(site == site_name)
    plots <- map2(pairs$group1, pairs$group2, function(g1, g2) {
        build_volcano(diffab_dir, site_name, g1, g2, label_n = label_n, title_site = FALSE) +
            theme(legend.position = "none")
    })
    shared_legend <- get_legend(
        build_volcano(diffab_dir, site_name, pairs$group1[1], pairs$group2[1],
                      label_n = label_n, title_site = FALSE) +
            theme(legend.position = "bottom")
    )
    n_rows <- ceiling(length(plots) / ncol)
    grid <- plot_grid(plotlist = plots, ncol = ncol)
    plot_grid(grid, shared_legend, ncol = 1, rel_heights = c(n_rows, 0.4))
}

## Setup
setwd(here::here())
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

## Differential abundance results must match whatever
## 9_differential_abundance_16s.R was run with (DIFF_AB_TEST_N for a fast
## test run, unset for the real one)
test_n <- suppressWarnings(as.integer(Sys.getenv("DIFF_AB_TEST_N", "")))
diffab_dir <- if (!is.na(test_n)) "results/differential_abundance_test" else "results/differential_abundance"
if (!is.na(test_n)) cat("TEST MODE: reading differential abundance results from", diffab_dir, "\n")

## ---- Panels A/B: pairwise DA taxa count heatmap, nose/throat ----
summary_pairs <- read_csv(file.path(diffab_dir, "summary_pairs_16s.csv"), show_col_types = FALSE)

## n_sig per pair: count qval < 0.05 in each pair's MaAsLin2 results (same
## threshold as the forest plots and the upset plot in 10_upset_diffabund_16s.R)
summary_pairs <- summary_pairs |>
    rowwise() |>
    mutate(n_sig = {
        f <- file.path(diffab_dir, "maaslin2",
                       paste0("maaslin2_ethnicity_results_16s_", site, "_",
                              pair_name(group1, group2), ".csv"))
        sum(read_csv(f, show_col_types = FALSE)$qval < 0.05, na.rm = TRUE)
    }) |>
    ungroup()

heatmap_panels <- list()
sites <- c("nose", "throat")

for (site_name in sites) {
    pair_counts <- summary_pairs |> filter(site == site_name)
    groups_order <- names(eth_colours)[names(eth_colours) %in%
                                        unique(c(pair_counts$group1, pair_counts$group2))]

    heatmap_panels[[site_name]] <- pairwise_sig_heatmap(
        pair_counts, groups_order,
        title = paste0("Pairwise DA: ", site_name),
        fill_limits = c(0, max(pair_counts$n_sig))
    )
}

## ---- Panel C: nose volcano highlight ----
## Dutch vs African Surinamese: best signal at this site (36/12); nose is
## underpowered everywhere else, so one panel represents it.
p_nose_1 <- build_volcano(diffab_dir, "nose", "Dutch", "African Surinamese")

## ---- Panels D/E/F: throat volcano highlights ----
## Dutch vs SAS: strongest signal (160 taxa q<0.25, 116 q<0.05). Dutch vs
## Ghanaian: second-strongest (106/72). Ghanaian vs Moroccan: strong signal
## (62/39) between two non-Dutch groups, showing migrant-migrant divergence
## is real, not just host-vs-migrant.
p_throat_1 <- build_volcano(diffab_dir, "throat", "Dutch", "South-Asian Surinamese")
p_throat_2 <- build_volcano(diffab_dir, "throat", "Dutch", "Ghanaian")
p_throat_3 <- build_volcano(diffab_dir, "throat", "Ghanaian", "Moroccan")

volcano_panels <- list(p_nose_1, p_throat_1, p_throat_2, p_throat_3)
shared_legend <- get_legend(volcano_panels[[1]] + theme(legend.position = "bottom"))
volcano_panels <- lapply(volcano_panels, function(p) p + theme(legend.position = "none"))

## ---- Assemble ----
grid <- plot_grid(
    heatmap_panels[["nose"]], heatmap_panels[["throat"]],
    volcano_panels[[1]], volcano_panels[[2]],
    volcano_panels[[3]], volcano_panels[[4]],
    ncol = 2, labels = c("A", "B", "C", "D", "E", "F"),
    rel_heights = c(0.8, 1.2, 1.2)
)

fig3 <- plot_grid(grid, shared_legend, ncol = 1, rel_heights = c(1, 0.04))

ggsave("results/figures/figure3.pdf", plot = fig3, width = 14, height = 16)
ggsave("results/figures/figure3.png", plot = fig3, width = 14, height = 16, dpi = 300)

cat("Saved results/figures/figure3.pdf and .png\n")

## ---- Supplement: every pairwise volcano plot, per site ----
## Not part of the lettered main figure - covers all 10 (nose) / 15 (throat)
## pairs, including the ones not highlighted above.
supp_nose <- assemble_site_volcanoes(diffab_dir, "nose", summary_pairs, ncol = 3)
n_rows_nose <- ceiling(sum(summary_pairs$site == "nose") / 3)
ggsave("results/figures/figure3_supp_nose_volcanoes.pdf", plot = supp_nose,
       width = 14, height = 4.5 * n_rows_nose + 1, limitsize = FALSE)
ggsave("results/figures/figure3_supp_nose_volcanoes.png", plot = supp_nose,
       width = 14, height = 4.5 * n_rows_nose + 1, dpi = 300, limitsize = FALSE)

supp_throat <- assemble_site_volcanoes(diffab_dir, "throat", summary_pairs, ncol = 3)
n_rows_throat <- ceiling(sum(summary_pairs$site == "throat") / 3)
ggsave("results/figures/figure3_supp_throat_volcanoes.pdf", plot = supp_throat,
       width = 14, height = 4.5 * n_rows_throat + 1, limitsize = FALSE)
ggsave("results/figures/figure3_supp_throat_volcanoes.png", plot = supp_throat,
       width = 14, height = 4.5 * n_rows_throat + 1, dpi = 300, limitsize = FALSE)

cat("Saved results/figures/figure3_supp_nose_volcanoes.pdf/.png and",
    "figure3_supp_throat_volcanoes.pdf/.png\n")

## Adjustment set per site, for the figure legend text (subtitles were
## dropped from the plots themselves)
for (site_name in sites) {
    cat(site_name, "- adjusted for:", site_confounder_label(diffab_dir, site_name), "\n")
}
