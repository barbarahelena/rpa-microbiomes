## Figure 1 assembly: 16S throat/nose, stratified by ethnicity
## Panels:
##   A/B - Shannon diversity (nose/throat)
##   C/D - Weighted UniFrac PCoA (nose/throat)
##   E/F - Betadisper, distance to centroid (nose/throat)
## Also writes a supplement (figure1_supp_braycurtis_pcoa.pdf/.png) with
## Bray-Curtis PCoA panels A/B for nose/throat side by side.
## Alpha diversity is recomputed directly from the rarefied phyloseq objects
## (cheap). The beta diversity panels are rebuilt from the .rds cache written
## by 7a_beta_diversity_16s_compute.R, so this script never repeats a
## permutation test - run 7a first (BETA_DIV_TEST_N=<n> for a fast test cache).

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(cowplot)
library(grid)
library(ggthemes)
library(ggpubr)

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

## Keep only ethnicity groups with more than n=50 samples (matches Table 1)
keep_groups <- function(ps, min_n = 50) {
    counts <- table(sample_data(ps)$EthnicityTotal)
    names(counts)[counts > min_n]
}

## PCoA ordination panel (points + 95% ellipses, centroids when >3 groups,
## PERMANOVA R2/p annotation) - shared by the main figure's Weighted UniFrac
## panels and the Bray-Curtis supplement.
build_pcoa_panel <- function(meta, pcoa, block, dist_name, site_name, eth_colours) {
    eig <- pcoa$values$Eigenvalues
    var_explained <- round(100 * eig / sum(eig), 1)

    ord_df <- data.frame(
        PCo1 = pcoa$vectors[, 1],
        PCo2 = pcoa$vectors[, 2],
        EthnicityTotal = meta$EthnicityTotal
    )

    permanova_label <- paste0(
        "R² = ", round(block$permanova_eth$R2[1], 3),
        ", p = ", format.pval(block$permanova_eth[["Pr(>F)"]][1], digits = 2, eps = 0.001)
    )

    n_groups <- nlevels(droplevels(meta$EthnicityTotal))

    p <- ggplot(ord_df, aes(x = PCo1, y = PCo2, colour = EthnicityTotal)) +
        geom_point(alpha = 0.5, size = 1) +
        stat_ellipse(level = 0.95, linewidth = 0.8)

    if (n_groups > 3) {
        centroids <- ord_df |>
            group_by(EthnicityTotal) |>
            summarise(PCo1 = mean(PCo1), PCo2 = mean(PCo2), .groups = "drop")
        p <- p +
            geom_point(data = centroids,
                       aes(x = PCo1, y = PCo2, fill = EthnicityTotal),
                       shape = 21, colour = "black", size = 4, stroke = 0.8) +
            scale_fill_manual(values = eth_colours, guide = "none")
    }

    p <- p +
        scale_colour_manual(values = eth_colours) +
        annotate("text", x = Inf, y = Inf, label = permanova_label,
                 hjust = 1.05, vjust = 1.5, size = 3.2) +
        labs(x = paste0("PCo1 (", var_explained[1], "%)"),
             y = paste0("PCo2 (", var_explained[2], "%)"),
             colour = "Ethnicity",
             title = paste0(dist_name, " - ", site_name)) +
        guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_Publication()

    list(plot = p, n_groups = n_groups)
}

## Setup
setwd(here::here())
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

## Beta diversity cache must match whatever 7a_beta_diversity_16s_compute.R
## was run with (BETA_DIV_TEST_N for a fast test cache, unset for the real one)
test_n <- suppressWarnings(as.integer(Sys.getenv("BETA_DIV_TEST_N", "")))
beta_outdir <- if (!is.na(test_n)) "results/beta_diversity_test" else "results/beta_diversity"
if (!is.na(test_n)) cat("TEST MODE: reading beta diversity cache from", beta_outdir, "\n")

## Ethnicity colours (same validated palette as scripts 6-8)
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

sites <- c("nose", "throat")

## ---- Panels A/B: Shannon diversity boxplots ----
shannon_panels <- list()

for (site_name in sites) {
    ps <- readRDS(paste0("data/processed/ps_", site_name, "_rarefied.RDS"))
    ps <- subset_samples(ps, EthnicityTotal %in% keep_groups(ps))

    alpha_df <- estimate_richness(ps, measures = "Shannon") |>
        rownames_to_column("sample_id")

    meta <- sample_data(ps) |>
        as("data.frame") |>
        rownames_to_column("sample_id") |>
        select(sample_id, EthnicityTotal) |>
        mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)))

    alpha_df <- alpha_df |>
        left_join(meta, by = "sample_id") |>
        filter(!is.na(EthnicityTotal))

    ## Kruskal-Wallis omnibus + pairwise Wilcoxon (BH-adjusted), same
    ## convention as 6_alpha_diversity_16s_ethnicity.R. Only p.adj < 0.05
    ## pairs get a bracket - up to 10 (nose) or 15 (throat) possible pairs
    ## would be unreadable if all were drawn.
    pw <- pairwise.wilcox.test(alpha_df$Shannon, alpha_df$EthnicityTotal, p.adjust.method = "BH")
    max_val <- max(alpha_df$Shannon)
    min_val <- min(alpha_df$Shannon)
    step <- (max_val - min_val) * 0.06

    sig_pairs <- as.data.frame(as.table(pw$p.value)) |>
        filter(!is.na(Freq)) |>
        rename(group1 = Var1, group2 = Var2, p.adj = Freq) |>
        filter(p.adj < 0.05) |>
        arrange(p.adj) |>
        mutate(
            group1 = as.character(group1),
            group2 = as.character(group2),
            y.position = max_val + step * row_number(),
            p.adj.label = case_when(
                p.adj < 0.0001 ~ "****",
                p.adj < 0.001  ~ "***",
                p.adj < 0.01   ~ "**",
                TRUE           ~ "*"
            )
        )

    p_shannon <- ggplot(alpha_df, aes(x = EthnicityTotal, y = Shannon,
                                       fill = EthnicityTotal)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
        stat_compare_means(method = "kruskal.test", label = "p.format") +
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = "Shannon index",
             title = paste0("Shannon diversity - ", site_name)) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        theme_Publication() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 35, hjust = 1))

    if (nrow(sig_pairs) > 0) {
        p_shannon <- p_shannon +
            ggpubr::stat_pvalue_manual(sig_pairs, label = "p.adj.label",
                                        xmin = "group1", xmax = "group2",
                                        y.position = "y.position",
                                        tip.length = 0, size = 3,
                                        color = "grey30")
    }

    shannon_panels[[site_name]] <- p_shannon
}

## ---- Panels C/D and E/F: PCoA (weighted UniFrac) and betadisper ----
pcoa_panels <- list()
betadisp_panels <- list()
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
    meta  <- cache$meta
    pcoa  <- cache$pcoas[["Weighted UniFrac"]]
    block <- cache$blocks[["Weighted UniFrac"]]

    ## PCoA
    res <- build_pcoa_panel(meta, pcoa, block, "Weighted UniFrac", site_name, eth_colours)
    p_pcoa <- res$plot
    n_groups <- res$n_groups

    ## The site with the most ethnicity groups gives the legend covering
    ## every colour used anywhere in the figure
    if (is.null(legend_source) || n_groups > legend_source$n_groups) {
        legend_source <- list(plot = p_pcoa, n_groups = n_groups)
    }

    pcoa_panels[[site_name]] <- p_pcoa + theme(legend.position = "none")

    ## Betadisper
    disp_df <- data.frame(
        Distance = block$betadisp$distances,
        EthnicityTotal = meta$EthnicityTotal
    )
    betadisp_panels[[site_name]] <- ggplot(disp_df, aes(x = EthnicityTotal, y = Distance,
                                                          fill = EthnicityTotal)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
        scale_fill_manual(values = eth_colours, guide = "none") +
        labs(x = NULL, y = "Distance to centroid",
             title = paste0("Betadisper - ", site_name),
             subtitle = paste0("Permutest p = ",
                               format.pval(block$betadisp_test$tab[["Pr(>F)"]][1], digits = 3))) +
        theme_Publication() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 35, hjust = 1))
}

shared_legend <- get_legend(legend_source$plot + theme(legend.position = "bottom"))

## ---- Assemble ----
grid <- plot_grid(
    shannon_panels[["nose"]], shannon_panels[["throat"]],
    pcoa_panels[["nose"]],    pcoa_panels[["throat"]],
    betadisp_panels[["nose"]], betadisp_panels[["throat"]],
    ncol = 2, labels = c("A", "B", "C", "D", "E", "F"),
    rel_heights = c(1.1, 1.2, 1.1)
)

fig1 <- plot_grid(grid, shared_legend, ncol = 1, rel_heights = c(1, 0.06))

ggsave("results/figures/figure1.pdf", plot = fig1, width = 10, height = 14)
ggsave("results/figures/figure1.png", plot = fig1, width = 10, height = 14, dpi = 300)

cat("Saved results/figures/figure1.pdf and .png\n")

## ---- Supplement: Bray-Curtis PCoA (nose/throat side by side) ----
## Reuses the cache objects already loaded above - no recompute.
bc_panels <- list()
bc_legend_source <- NULL

for (site_name in sites) {
    cache <- caches[[site_name]]
    meta  <- cache$meta
    pcoa  <- cache$pcoas[["Bray-Curtis"]]
    block <- cache$blocks[["Bray-Curtis"]]

    res <- build_pcoa_panel(meta, pcoa, block, "Bray-Curtis", site_name, eth_colours)

    if (is.null(bc_legend_source) || res$n_groups > bc_legend_source$n_groups) {
        bc_legend_source <- list(plot = res$plot, n_groups = res$n_groups)
    }

    bc_panels[[site_name]] <- res$plot + theme(legend.position = "none")
}

bc_shared_legend <- get_legend(bc_legend_source$plot + theme(legend.position = "bottom"))

bc_grid <- plot_grid(
    bc_panels[["nose"]], bc_panels[["throat"]],
    ncol = 2, labels = c("A", "B")
)

fig1_supp_bc <- plot_grid(bc_grid, bc_shared_legend, ncol = 1, rel_heights = c(1, 0.12))

ggsave("results/figures/figure1_supp_braycurtis_pcoa.pdf", plot = fig1_supp_bc, width = 10, height = 5.5)
ggsave("results/figures/figure1_supp_braycurtis_pcoa.png", plot = fig1_supp_bc, width = 10, height = 5.5, dpi = 300)

cat("Saved results/figures/figure1_supp_braycurtis_pcoa.pdf and .png\n")
