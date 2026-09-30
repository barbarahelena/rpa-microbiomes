## Compute participant air-pollution summaries and tests.

## Libraries
library(here)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(phyloseq)

## Setup
setwd(here::here())
dir.create("results/airpollution", recursive = TRUE, showWarnings = FALSE)

## Keep only ethnicity groups with more than n=50 samples (matches Table 1)
drop_small_groups <- function(data, min_n = 50) {
    data |>
        add_count(EthnicityTotal, name = "n_group") |>
        filter(n_group > min_n) |>
        select(-n_group) |>
        droplevels()
}

## Data
meta <- readRDS("data/processed/HELIUSmetadata_clean.RDS") |>
    filter(!is.na(PM25_mean)) |>
    drop_small_groups()
cat("Participants with air pollution data:", nrow(meta), "\n")

pollutants <- c(
    PM10_mean = "PM10 (µg/m³, 2013-2015 mean)",
    PM25_mean = "PM2.5 (µg/m³, 2013-2015 mean)",
    NO2_mean  = "NO2 (µg/m³, 2014-2015 mean)",
    EC_mean   = "Soot/EC (µg/m³, 2013-2015 mean)"
)

## Pairwise Wilcoxon (BH-adjusted) per pollutant - kept in full in the CSV;
## only p.adj < 0.05 pairs get bracket annotations on the boxplots (with up
## to 8 ethnicity groups / 28 pairs, showing every pair would be unreadable).
pairwise_results <- lapply(names(pollutants), function(var) {
    pw <- pairwise.wilcox.test(meta[[var]], meta$EthnicityTotal, p.adjust.method = "BH")
    as.data.frame(as.table(pw$p.value)) |>
        filter(!is.na(Freq)) |>
        dplyr::rename(group1 = Var1, group2 = Var2, p.adj = Freq) |>
        mutate(pollutant = var, .before = 1)
}) |> bind_rows()
## Summary table: mean/median/range per pollutant, overall and by ethnicity
summary_overall <- meta |>
    summarise(across(all_of(names(pollutants)),
                      list(mean = ~mean(.x), median = ~median(.x),
                           min = ~min(.x), max = ~max(.x)))) |>
    mutate(EthnicityTotal = "All participants", .before = 1)

summary_by_eth <- meta |>
    group_by(EthnicityTotal) |>
    summarise(across(all_of(names(pollutants)),
                      list(mean = ~mean(.x), median = ~median(.x),
                           min = ~min(.x), max = ~max(.x)))) |>
    mutate(EthnicityTotal = as.character(EthnicityTotal))

## Tests previously evaluated by the plot layer; keep its formatting convention.
kruskal_annotations <- lapply(names(pollutants), function(var) {
    ggpubr::compare_means(reformulate("EthnicityTotal", response = var),
                         data = meta, method = "kruskal.test")
}) |> setNames(names(pollutants))

## Figure 1 orders groups by median before its pairwise tests. Preserve that
## ordering here because it determines pair orientation and bracket ordering.
figure_tests <- lapply(names(pollutants), function(var) {
    df <- meta |> select(EthnicityTotal, value = all_of(var))
    ordered_levels <- df |>
        group_by(EthnicityTotal) |>
        summarise(m = median(value), .groups = "drop") |>
        arrange(m) |>
        pull(EthnicityTotal) |>
        as.character()
    df <- df |> mutate(EthnicityTotal = factor(EthnicityTotal, levels = ordered_levels))
    list(kw_p = kruskal.test(value ~ EthnicityTotal, data = df)$p.value,
         pw = pairwise.wilcox.test(df$value, df$EthnicityTotal, p.adjust.method = "BH"))
}) |> setNames(names(pollutants))

dir.create("results/airpollution/cache", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(meta = meta, pollutants = pollutants, pairwise_results = pairwise_results,
        summary_overall = summary_overall, summary_by_eth = summary_by_eth,
        kruskal_annotations = kruskal_annotations, figure_tests = figure_tests),
        "results/airpollution/cache/participant_exposure.rds")
