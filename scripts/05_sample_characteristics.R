## Sampling seasonality and sequencing-batch composition.

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

## ---- Seasonality of swab collection by ethnicity ----
## Day-of-year for the 1st of each month (non-leap reference year), used as
## x-axis gridlines/labels so the plot reads by calendar month
month_starts <- yday(as.Date(paste0("2001-", 1:12, "-01")))

## Unrarefied, QC'd phyloseq objects (post decontam/dedup, pre rarefaction) -
## rarefaction-driven sample dropout isn't relevant to a sampling-date check
sites <- list(
    throat = readRDS("data/processed/ps_throat.RDS"),
    nose   = readRDS("data/processed/ps_nose.RDS")
)

for (site_name in names(sites)) {
    ps <- sites[[site_name]]

    ## Filter to ethnicity groups with N > 50 in this site (matches above)
    keep <- table(sample_data(ps)$EthnicityTotal)
    ps <- subset_samples(ps, EthnicityTotal %in% names(keep)[keep > 50])

    date_df <- sample_data(ps) |>
        as("data.frame") |>
        filter(!is.na(Collection_Date)) |>
        mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)),
               yday = yday(Collection_Date))

    dir.create("results/sample_metadata/cache", recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(date_df = date_df),
            paste0("results/sample_metadata/cache/seasonality_16s_", site_name, ".rds"))

    n_missing <- nsamples(ps) - nrow(date_df)
    cat(site_name, ": dropped", n_missing, "sample(s) with missing collection date\n")

    ## Monthly visit counts by ethnicity (reference table)
    month_counts <- date_df |>
        mutate(month = month(Collection_Date, label = TRUE)) |>
        count(EthnicityTotal, month, name = "n") |>
        complete(EthnicityTotal, month, fill = list(n = 0))
    write_csv(month_counts,
              paste0("results/sample_metadata/season_by_ethnicity_16s_", site_name, "_counts.csv"))

    ## Density plot: day-of-year sampling distribution by ethnicity
    p <- ggplot(date_df, aes(x = yday, colour = EthnicityTotal, fill = EthnicityTotal)) +
        geom_density(alpha = 0.15, linewidth = 0.8) +
        scale_colour_manual(values = eth_colours, name = "Ethnicity") +
        scale_fill_manual(values = eth_colours, name = "Ethnicity") +
        scale_x_continuous(breaks = month_starts, labels = month.abb,
                            limits = c(1, 366), expand = c(0, 0)) +
        labs(title = paste0("Sampling seasonality by ethnicity - 16S ", site_name),
             x = "Collection month", y = "Density") +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(paste0("results/sample_metadata/season_by_ethnicity_16s_", site_name, ".pdf"),
           p, width = 8, height = 5)
    ggsave(paste0("results/sample_metadata/season_by_ethnicity_16s_", site_name, ".png"),
           p, width = 8, height = 5, dpi = 300)
}

## ---- Batch composition by ethnicity ----
## Distribution of ethnicity groups across sequencing runs (SeqBatch) -
## context for the batch covariate used in the alpha/beta diversity and
## differential abundance analyses. DNAIsoBatch (DNA isolation
## date) was also tested as a candidate batch covariate but explains largely
## overlapping beta-diversity variance and isn't used (see
## scripts/10_beta_diversity_16s_ethnicity_compute.R), so it has no distribution plot here.
batch_vars <- c(
    SeqBatch = "Sequencing batch"
)

for (site_name in names(sites)) {
    ps <- sites[[site_name]]

    ## Filter to ethnicity groups with N > 50 in this site (matches above)
    keep <- table(sample_data(ps)$EthnicityTotal)
    ps <- subset_samples(ps, EthnicityTotal %in% names(keep)[keep > 50])
    site_meta <- sample_data(ps) |>
        as("data.frame") |>
        mutate(EthnicityTotal = droplevels(factor(EthnicityTotal)))

    for (batch_var in names(batch_vars)) {
        batch_df <- site_meta |>
            filter(!is.na(.data[[batch_var]])) |>
            rename(Batch = all_of(batch_var))

        n_missing <- nrow(site_meta) - nrow(batch_df)
        cat(site_name, "-", batch_var, ": dropped", n_missing,
            "sample(s) with missing batch\n")

        ## Batch labels (isolation dates / sequential run IDs) sort
        ## chronologically as plain strings
        batch_df$Batch <- factor(batch_df$Batch, levels = sort(unique(batch_df$Batch)))

        ## Ethnicity counts per batch (reference table)
        batch_counts <- batch_df |>
            count(Batch, EthnicityTotal, name = "n") |>
            complete(Batch, EthnicityTotal, fill = list(n = 0))
        write_csv(batch_counts,
                  paste0("results/sample_metadata/", tolower(batch_var), "_by_ethnicity_16s_",
                         site_name, "_counts.csv"))

        ## Stacked proportional bar: ethnicity composition per batch
        p <- ggplot(batch_df, aes(x = Batch, fill = EthnicityTotal)) +
            geom_bar(position = "fill") +
            scale_fill_manual(values = eth_colours, name = "Ethnicity") +
            scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
            labs(title = paste0(batch_vars[[batch_var]], " composition by ethnicity - 16S ",
                                 site_name),
                 x = batch_vars[[batch_var]], y = "Proportion of samples") +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

        ggsave(paste0("results/sample_metadata/", tolower(batch_var), "_by_ethnicity_16s_",
                      site_name, ".pdf"),
               p, width = 10, height = 5)
        ggsave(paste0("results/sample_metadata/", tolower(batch_var), "_by_ethnicity_16s_",
                      site_name, ".png"),
               p, width = 10, height = 5, dpi = 300)
    }
}
