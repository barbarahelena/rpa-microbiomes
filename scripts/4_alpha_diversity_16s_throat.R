## Alpha diversity analysis: 16S throat microbiome
## Descriptive only — no ethnicity stratification (no metadata linkage available)

## Libraries
library(here)
library(tidyverse)
library(phyloseq)

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
dir.create("results/alpha_diversity", recursive = TRUE, showWarnings = FALSE)

## Load rarefied phyloseq object
ps_throat <- readRDS("data/processed/ps_throat_rarefied.RDS")

## Compute alpha diversity metrics
alpha_df <- estimate_richness(ps_throat, measures = c("Observed", "Shannon", "Simpson")) |>
    rownames_to_column("sample_id")

## Pivot to long format for summaries and faceted plots
alpha_long <- alpha_df |>
    pivot_longer(cols = c(Observed, Shannon, Simpson),
                 names_to = "metric", values_to = "value")

## Summary statistics table
summary_table <- alpha_long |>
    group_by(metric) |>
    summarise(
        n      = n(),
        mean   = mean(value),
        sd     = sd(value),
        median = median(value),
        q25    = quantile(value, 0.25),
        q75    = quantile(value, 0.75),
        min    = min(value),
        max    = max(value),
        .groups = "drop"
    )
write_csv(summary_table, "results/alpha_diversity/alpha_diversity_summary_16s_throat.csv")

## Figure 1: Violin + boxplot (primary figure)
n_samples <- nsamples(ps_throat)
ggplot(alpha_long, aes(x = metric, y = value)) +
    geom_violin(fill = "#A6CEE3", alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = 21, outlier.size = 0.8) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = "Value",
         title = paste0("Alpha diversity - 16S throat microbiome (n = ", n_samples, ")"),
         caption = "Rarefied to 7,500 reads/sample") +
    theme_Publication() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
ggsave("results/alpha_diversity/alpha_diversity_violin_16s_throat.pdf",
       width = 9, height = 5)

## Figure 2: Histograms (supplementary)
ggplot(alpha_long, aes(x = value)) +
    geom_histogram(fill = "#A6CEE3", colour = "black", linewidth = 0.2, bins = 40) +
    facet_wrap(~ metric, scales = "free", nrow = 1) +
    labs(x = "Value", y = "Count",
         title = "Alpha diversity distributions - 16S throat",
         caption = "Rarefied to 7,500 reads/sample") +
    theme_Publication()
ggsave("results/alpha_diversity/alpha_diversity_hist_16s_throat.pdf",
       width = 9, height = 4)

## Figure 3: Pairwise scatter (supplementary)
ggplot(alpha_df, aes(x = Observed, y = Shannon, colour = Simpson)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_colour_viridis_c(option = "plasma") +
    labs(x = "Observed species", y = "Shannon index", colour = "Simpson",
         title = "Alpha diversity metric relationships - 16S throat",
         caption = "Rarefied to 7,500 reads/sample") +
    theme_Publication() +
    theme(legend.position = "right")
ggsave("results/alpha_diversity/alpha_diversity_scatter_16s_throat.pdf",
       width = 7, height = 5)
