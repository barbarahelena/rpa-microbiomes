## Report participant air-pollution exposure from saved tests.

## Libraries
library(here)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(phyloseq)

## Functions
source(here::here("scripts", "lib", "plot_style.R"))

## Setup
setwd(here::here())
dir.create("results/airpollution", recursive = TRUE, showWarnings = FALSE)
dir.create("results/sample_metadata", recursive = TRUE, showWarnings = FALSE)

## Shared ethnicity colours
eth_colours <- ethnicity_colours

cache_path <- "results/airpollution/cache/participant_exposure.rds"
if (!file.exists(cache_path)) stop("Missing cache: ", cache_path, "; run pixi run airpollution-participants-compute first.")
cache <- readRDS(cache_path)
meta <- cache$meta
pollutants <- cache$pollutants
pairwise_results <- cache$pairwise_results
summary_overall <- cache$summary_overall
summary_by_eth <- cache$summary_by_eth
kruskal_annotations <- cache$kruskal_annotations

write_csv(pairwise_results, "results/airpollution/participants_pairwise_wilcoxon.csv")

## One histogram + one ethnicity boxplot per pollutant
## Vertical headroom above the boxes (Kruskal-Wallis p-value + stacked
## pairwise brackets) is sized to the number of significant pairs found, via
## expand(mult = ...) on the y-axis, so labels never get clipped at the top.
plots <- lapply(names(pollutants), function(var) {
    hist <- ggplot(meta, aes(x = .data[[var]])) +
        geom_histogram(fill = "steelblue", bins = 40) +
        labs(x = pollutants[[var]], y = "Participants", title = pollutants[[var]]) +
        theme_Publication()

    max_val <- max(meta[[var]])
    min_val <- min(meta[[var]])
    step <- (max_val - min_val) * 0.12

    ## Brackets stack from just above the boxes upward, most-significant pair
    ## first; the Kruskal-Wallis omnibus label sits above all of them so it
    ## never collides with a bracket line. Capped at the 6 most significant
    ## pairs - beyond that the brackets overlap and become unreadable (full
    ## pairwise results, capped or not, are always in the CSV).
    sig_pairs <- pairwise_results |>
        filter(pollutant == var, p.adj < 0.05) |>
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
    n_brackets <- nrow(sig_pairs)
    kruskal_y <- max_val + step * (n_brackets + 1.5)

    ## Bracket y-positions are real data points that ggplot trains the y-scale
    ## to, so the panel already extends to cover them - a small fixed expand
    ## on top of that (not scaled by n_brackets) is enough to keep the top
    ## label clear of the edge without squeezing the boxes down to a sliver.
    box <- ggplot(meta, aes(x = EthnicityTotal, y = .data[[var]], fill = EthnicityTotal)) +
        geom_boxplot(outlier.size = 0.8) +
        stat_cached_pvalue(paste0("p = ", kruskal_annotations[[var]]$p.format),
                           label.y = kruskal_y) +
        scale_fill_manual(values = eth_colours, guide = "none") +
        labs(x = NULL, y = pollutants[[var]]) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 40, hjust = 1))

    if (n_brackets > 0) {
        box <- box +
            stat_pvalue_manual(sig_pairs, label = "p.adj.label",
                                xmin = "group1", xmax = "group2",
                                y.position = "y.position",
                                tip.length = 0.01, size = 3)
    }

    ggarrange(hist, box, nrow = 1)
})

for (i in seq_along(pollutants)) {
    ggsave(
        paste0("results/airpollution/participants_", names(pollutants)[i], ".pdf"),
        plots[[i]], width = 10, height = 4.5
    )
}

ggarrange(plotlist = plots, ncol = 1)
ggsave("results/airpollution/participants_all_pollutants.pdf", width = 10, height = 16)

bind_rows(summary_overall, summary_by_eth) |>
    write_csv("results/airpollution/participants_summary.csv")

