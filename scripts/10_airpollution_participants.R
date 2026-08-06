## Air pollution exposure among HELIUS participants (this cohort only)
## Distribution of individual RIVM/ALO exposure estimates (2013-2015 mean,
## see 1a_datacleaning_helius.R for the Heliusnr-based linkage), overall and
## by ethnicity. Participant addresses/coordinates are not available (GECCO
## delivers pre-extracted per-participant values only, not addresses), so
## this is a distribution plot rather than a map - the city-wide PC6 map is
## in 11_airpollution_amsterdam.R.

## Libraries
library(here)
library(tidyverse)
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

## Setup
setwd(here::here())
dir.create("results/airpollution", recursive = TRUE, showWarnings = FALSE)

## Ethnicity colours (same validated palette as scripts 4, 7 and 8)
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

## One histogram + one ethnicity boxplot per pollutant
plots <- lapply(names(pollutants), function(var) {
    hist <- ggplot(meta, aes(x = .data[[var]])) +
        geom_histogram(fill = "steelblue", bins = 40) +
        labs(x = pollutants[[var]], y = "Participants", title = pollutants[[var]]) +
        theme_Publication()

    box <- ggplot(meta, aes(x = EthnicityTotal, y = .data[[var]], fill = EthnicityTotal)) +
        geom_boxplot(outlier.size = 0.8) +
        scale_fill_manual(values = eth_colours, guide = "none") +
        labs(x = NULL, y = pollutants[[var]]) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 40, hjust = 1))

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

bind_rows(summary_overall, summary_by_eth) |>
    write_csv("results/airpollution/participants_summary.csv")
