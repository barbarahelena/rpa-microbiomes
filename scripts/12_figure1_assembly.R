## Figure 1 assembly: 16S throat/nose, stratified by ethnicity
## Panels:
##   A/B/C/D - Amsterdam air pollution maps (PM10, PM2.5, NO2, EC; small, static)
##   E/F/G/H - Participant-level air pollution exposure by ethnicity (PM10,
##             PM2.5, NO2, EC boxplots; each x-axis ordered by group median)
##   I/J/K/L - Sampling seasonality density (I/J, nose above throat) and
##             Shannon diversity (K/L, nose above throat) as two columns
##             side by side, density wider and Shannon narrower to match
##             their content
##   M/N     - Weighted UniFrac PCoA (nose/throat)
## Also writes:
##   - figure1_supp_braycurtis_pcoa.pdf/.png: Bray-Curtis PCoA panels A/B
##     for nose/throat side by side.
##   - figure1_supp_betadisper.pdf/.png: betadisper (distance to centroid)
##     panels A/B for nose/throat side by side.
## Alpha diversity, the density panels, and the map are recomputed directly
## from the underlying phyloseq/geometry objects (cheap). The beta diversity
## panels (PCoA, betadisper) are rebuilt from the .rds cache written by
## 7a_beta_diversity_16s_compute.R, so this script never repeats a
## permutation test - run 7a first (BETA_DIV_TEST_N=<n> for a fast test cache).
## The map panel is rebuilt from the .rds geometry cache written by
## 4_airpollution_amsterdam.R - run that first too.

## Libraries
library(here)
library(tidyverse)
library(phyloseq)
library(cowplot)
library(grid)
library(ggthemes)
library(ggpubr)
library(sf)

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

## Same threshold, for a plain data.frame (participant-level metadata) rather
## than a phyloseq object
keep_groups_df <- function(df, min_n = 50) {
    counts <- table(df$EthnicityTotal)
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

    ## PERMANOVA R2/p goes in the subtitle (above the plot area) rather
    ## than as an in-plot annotation, so it doesn't sit on top of points
    p <- p +
        scale_colour_manual(values = eth_colours) +
        labs(x = paste0("PCo1 (", var_explained[1], "%)"),
             y = paste0("PCo2 (", var_explained[2], "%)"),
             colour = "Ethnicity",
             title = paste0(dist_name, " - ", site_name),
             subtitle = permanova_label) +
        guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_Publication() +
        theme(plot.subtitle = element_text(size = rel(0.75), hjust = 0.5))

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

## ---- Panels A/B/C/D: Amsterdam air pollution maps ----
geo_cache_path <- "results/airpollution/amsterdam_pc6_geo.rds"
if (!file.exists(geo_cache_path)) {
    stop("No Amsterdam PC6 geometry cache at ", geo_cache_path,
         " - run scripts/4_airpollution_amsterdam.R first.")
}
geo <- readRDS(geo_cache_path)

## coord_sf(expand = FALSE) plus zero plot margins keep each map filling its
## panel instead of floating in whitespace (see 4_airpollution_amsterdam.R
## for why the source geometries can otherwise produce a wildly oversized
## bounding box)
map_specs <- list(
    list(col = "pm10_avg_2013_2015", title = "PM10 (2013-2015 mean)",    legend = "PM10\n(µg/m³)"),
    list(col = "pm25_avg_2013_2015", title = "PM2.5 (2013-2015 mean)",   legend = "PM2.5\n(µg/m³)"),
    list(col = "no2_avg_2014_2015",  title = "NO2 (2014-2015 mean)",     legend = "NO2\n(µg/m³)"),
    list(col = "ec_avg_2013_2015",   title = "Soot/EC (2013-2015 mean)", legend = "EC\n(µg/m³)")
)

build_map_panel <- function(spec, geo) {
    ggplot(geo$pc6) +
        geom_sf(aes(fill = .data[[spec$col]]), colour = NA) +
        geom_sf(data = geo$boundary, fill = NA, colour = "black", linewidth = 0.3) +
        coord_sf(expand = FALSE) +
        scale_fill_viridis_c(name = spec$legend, option = "magma", direction = -1,
                              na.value = "grey85",
                              guide = guide_colorbar(barwidth = unit(0.2, "cm"), barheight = unit(1.5, "cm"))) +
        labs(title = spec$title) +
        theme_void(base_size = 10) +
        theme(plot.title = element_text(face = "bold", size = rel(0.9), hjust = 0.5),
              plot.margin = margin(2, 2, 2, 2),
              legend.text = element_text(size = rel(0.6)),
              legend.title = element_text(size = rel(0.7), face = "bold"))
}

map_panels <- lapply(map_specs, build_map_panel, geo = geo)

## ---- Panels E/F/G/H: participant-level air pollution exposure by ethnicity ----
## Unlike the city-wide PC6 maps above, these are per-participant RIVM/ALO
## estimates linked to each HELIUS participant (see 3_airpollution_participants.R,
## which this reuses the exact data/threshold/test convention from). Each
## panel's x-axis is reordered by that pollutant's group mean (low to high)
## rather than a fixed ethnicity order, so the gradient reads directly
## without cross-referencing the legend.
pollution_meta <- readRDS("data/processed/HELIUSmetadata_clean.RDS") |>
    filter(!is.na(PM25_mean))
pollution_meta <- pollution_meta |>
    filter(EthnicityTotal %in% keep_groups_df(pollution_meta)) |>
    mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)))

pollution_specs <- list(
    list(col = "PM10_mean", title = "PM10 (2013-2015 mean)",    y = "PM10 (µg/m³)"),
    list(col = "PM25_mean", title = "PM2.5 (2013-2015 mean)",   y = "PM2.5 (µg/m³)"),
    list(col = "NO2_mean",  title = "NO2 (2014-2015 mean)",     y = "NO2 (µg/m³)"),
    list(col = "EC_mean",   title = "Soot/EC (2013-2015 mean)", y = "Soot/EC (µg/m³)")
)

build_pollution_panel <- function(spec, meta, eth_colours) {
    df <- meta |> select(EthnicityTotal, value = all_of(spec$col))

    ## Order ethnicity groups by median exposure, low to high
    ordered_levels <- df |>
        group_by(EthnicityTotal) |>
        summarise(m = median(value), .groups = "drop") |>
        arrange(m) |>
        pull(EthnicityTotal) |>
        as.character()
    df <- df |> mutate(EthnicityTotal = factor(EthnicityTotal, levels = ordered_levels))

    ## Kruskal-Wallis omnibus + pairwise Wilcoxon (BH-adjusted), same
    ## convention as the Shannon panels above - only the 6 most significant
    ## pairs get a bracket, full results aren't written here since this
    ## duplicates 3_airpollution_participants.R's own CSV output. The
    ## omnibus p-value is placed in the subtitle (above the plot area)
    ## rather than drawn onto the panel with stat_compare_means, so it
    ## doesn't compete for space with the boxes/brackets.
    kw_p <- kruskal.test(value ~ EthnicityTotal, data = df)$p.value
    pw <- pairwise.wilcox.test(df$value, df$EthnicityTotal, p.adjust.method = "BH")
    max_val <- max(df$value)
    min_val <- min(df$value)
    step <- (max_val - min_val) * 0.08

    sig_pairs <- as.data.frame(as.table(pw$p.value)) |>
        filter(!is.na(Freq)) |>
        dplyr::rename(group1 = Var1, group2 = Var2, p.adj = Freq) |>
        filter(p.adj < 0.05) |>
        arrange(p.adj) |>
        slice_head(n = 6) |>
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

    p <- ggplot(df, aes(x = EthnicityTotal, y = value, fill = EthnicityTotal)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 0.8, alpha = 0.7) +
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = spec$y, title = spec$title,
             subtitle = paste0("Kruskal-Wallis p ", format.pval(kw_p, digits = 2, eps = 0.0001))) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        theme_Publication() +
        theme(legend.position = "none",
              plot.subtitle = element_text(size = rel(0.75), hjust = 0.5),
              axis.text.x = element_text(angle = 35, hjust = 1))

    if (nrow(sig_pairs) > 0) {
        p <- p +
            ggpubr::stat_pvalue_manual(sig_pairs, label = "p.adj.label",
                                        xmin = "group1", xmax = "group2",
                                        y.position = "y.position",
                                        tip.length = 0, size = 3, color = "grey30")
    }
    p
}

pollution_panels <- lapply(pollution_specs, build_pollution_panel,
                            meta = pollution_meta, eth_colours = eth_colours)

## ---- Panels I/J: sampling seasonality density plots (nose/throat) ----
## Day-of-year for the 1st of each month (non-leap reference year), used as
## x-axis gridlines/labels so the plot reads by calendar month
month_starts <- yday(as.Date(paste0("2001-", 1:12, "-01")))

build_density_panel <- function(site_name, eth_colours) {
    ## Unrarefied, QC'd phyloseq object (post decontam/dedup, pre
    ## rarefaction) - rarefaction-driven sample dropout isn't relevant to a
    ## sampling-date check
    ps <- readRDS(paste0("data/processed/ps_", site_name, ".RDS"))
    ## subset_samples()'s non-standard evaluation of its subset expression
    ## doesn't reliably resolve local variables when called from inside a
    ## function invoked via lapply (unlike the top-level loops below) -
    ## prune_samples() takes a plain logical vector instead, sidestepping
    ## the NSE lookup entirely
    keep <- keep_groups(ps)
    ps <- prune_samples(sample_data(ps)$EthnicityTotal %in% keep, ps)

    date_df <- sample_data(ps) |>
        as("data.frame") |>
        filter(!is.na(Collection_Date)) |>
        mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)),
               yday = yday(Collection_Date))

    ggplot(date_df, aes(x = yday, colour = EthnicityTotal, fill = EthnicityTotal)) +
        geom_density(alpha = 0.15, linewidth = 0.8) +
        scale_colour_manual(values = eth_colours) +
        scale_fill_manual(values = eth_colours) +
        scale_x_continuous(breaks = month_starts, labels = month.abb,
                            limits = c(1, 366), expand = c(0, 0)) +
        labs(x = "Collection month", y = "Density",
             title = paste0("Sampling season - ", site_name)) +
        theme_Publication() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))
}

density_panels <- lapply(sites, build_density_panel, eth_colours = eth_colours) |>
    setNames(sites)

## ---- Panels K/L: Shannon diversity boxplots ----
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
    ## would be unreadable if all were drawn. The omnibus p-value is placed
    ## in the subtitle (above the plot area) rather than drawn onto the
    ## panel with stat_compare_means, so it doesn't compete for space with
    ## the boxes/brackets.
    kw_p <- kruskal.test(Shannon ~ EthnicityTotal, data = alpha_df)$p.value
    pw <- pairwise.wilcox.test(alpha_df$Shannon, alpha_df$EthnicityTotal, p.adjust.method = "BH")
    max_val <- max(alpha_df$Shannon)
    min_val <- min(alpha_df$Shannon)
    step <- (max_val - min_val) * 0.06

    sig_pairs <- as.data.frame(as.table(pw$p.value)) |>
        filter(!is.na(Freq)) |>
        dplyr::rename(group1 = Var1, group2 = Var2, p.adj = Freq) |>
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
        scale_fill_manual(values = eth_colours) +
        labs(x = NULL, y = "Shannon index",
             title = paste0("Shannon diversity - ", site_name),
             subtitle = paste0("Kruskal-Wallis p ", format.pval(kw_p, digits = 2, eps = 0.0001))) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        theme_Publication() +
        theme(legend.position = "none",
              plot.subtitle = element_text(size = rel(0.75), hjust = 0.5),
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

## ---- Panels M/N (PCoA) and betadisper supplement panels A/B ----
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
row_maps <- plot_grid(
    plotlist = map_panels,
    ncol = 4, labels = c("A", "B", "C", "D")
)

row_pollution <- plot_grid(
    plotlist = pollution_panels,
    ncol = 4, labels = c("E", "F", "G", "H")
)

## Density (nose/throat stacked) and Shannon (nose/throat stacked) form two
## columns side by side, density wider and Shannon narrower to match the
## shape that suits each: density panels wider than tall, Shannon panels
## taller than wide.
col_density <- plot_grid(
    density_panels[["nose"]], density_panels[["throat"]],
    ncol = 1, labels = c("I", "J")
)

col_shannon <- plot_grid(
    shannon_panels[["nose"]], shannon_panels[["throat"]],
    ncol = 1, labels = c("K", "L")
)

row_density_shannon <- plot_grid(
    col_density, col_shannon,
    ncol = 2, rel_widths = c(1.3, 0.4)
)

row_pcoa <- plot_grid(
    pcoa_panels[["nose"]], pcoa_panels[["throat"]],
    ncol = 2, labels = c("M", "N")
)

grid <- plot_grid(row_maps, row_pollution, row_density_shannon, row_pcoa,
                   ncol = 1, rel_heights = c(0.7, 1.2, 1.9, 1.5))

fig1 <- plot_grid(grid, shared_legend, ncol = 1, rel_heights = c(1, 0.03))

ggsave("results/figures/figure1.pdf", plot = fig1, width = 15, height = 20)
ggsave("results/figures/figure1.png", plot = fig1, width = 15, height = 20, dpi = 300)

cat("Saved results/figures/figure1.pdf and .png\n")

## ---- Supplement: betadisper, distance to centroid (nose/throat) ----
betadisp_grid <- plot_grid(
    betadisp_panels[["nose"]], betadisp_panels[["throat"]],
    ncol = 2, labels = c("A", "B")
)

fig1_supp_betadisper <- plot_grid(betadisp_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.12))

ggsave("results/figures/figure1_supp_betadisper.pdf", plot = fig1_supp_betadisper, width = 10, height = 5.5)
ggsave("results/figures/figure1_supp_betadisper.png", plot = fig1_supp_betadisper, width = 10, height = 5.5, dpi = 300)

cat("Saved results/figures/figure1_supp_betadisper.pdf and .png\n")

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
