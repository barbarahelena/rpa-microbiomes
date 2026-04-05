## Descriptive composition plots: stacked bar plots of relative abundance
## 16S: overall composition (no ethnicity metadata available)
## Shotgun: split by ethnicity (Dutch vs South-Asian Surinamese)

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

# Helper: summarise mean relative abundance, keep top_n taxa, lump rest as "Other"
# df must have columns: sample_id, Abundance, and the tax_rank column
prep_abundance <- function(df, tax_rank, top_n, group_col = NULL) {
  df <- df |> mutate(taxon = replace_na(.data[[tax_rank]], "Unclassified"))

  # First: sum abundance per taxon per sample (collapse ASVs/species into taxon)
  group_vars <- c("sample_id", "taxon")
  if (!is.null(group_col)) group_vars <- c(group_vars, group_col)
  df_per_sample <- df |>
    group_by(across(all_of(group_vars))) |>
    summarise(Abundance = sum(Abundance), .groups = "drop")

  # Determine top taxa by mean abundance across all samples
  top_taxa <- df_per_sample |>
    group_by(taxon) |>
    summarise(mean_abund = mean(Abundance), .groups = "drop") |>
    arrange(desc(mean_abund)) |>
    slice_head(n = top_n) |>
    pull(taxon)

  df_per_sample <- df_per_sample |>
    mutate(taxon = if_else(taxon %in% top_taxa, taxon, "Other")) |>
    group_by(across(all_of(group_vars))) |>
    summarise(Abundance = sum(Abundance), .groups = "drop")

  # Mean across samples (per group if provided)
  if (!is.null(group_col)) {
    df_summary <- df_per_sample |>
      group_by(.data[[group_col]], taxon) |>
      summarise(mean_abund = mean(Abundance), .groups = "drop")
  } else {
    df_summary <- df_per_sample |>
      group_by(taxon) |>
      summarise(mean_abund = mean(Abundance), .groups = "drop")
  }

  # Order: largest taxa first, "Other" and "Unclassified" last
  tax_order <- df_summary |>
    group_by(taxon) |>
    summarise(total = sum(mean_abund), .groups = "drop") |>
    arrange(desc(total)) |>
    pull(taxon)
  tax_order <- c(setdiff(tax_order, c("Other", "Unclassified")),
                 intersect(tax_order, c("Unclassified", "Other")))
  df_summary <- df_summary |> mutate(taxon = factor(taxon, levels = tax_order))

  df_summary
}

# Set working directory to project root and create output directories if needed
setwd(here::here())
dir.create("results/composition", recursive = TRUE, showWarnings = FALSE)

# ---- 16S composition plots ----
ps_nose <- readRDS("data/processed/ps_nose_rarefied.RDS")
ps_throat <- readRDS("data/processed/ps_throat_rarefied.RDS")

process_16s <- function(ps, site_label) {
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  df <- psmelt(ps_rel) |> rename(sample_id = Sample)
  df$site <- site_label
  df
}

df_16s <- bind_rows(
  process_16s(ps_nose, "Nose"),
  process_16s(ps_throat, "Throat")
)

# Phylum-level plot
df_16s_phylum <- prep_abundance(df_16s, "Phylum", top_n = 10, group_col = "site")

n_taxa <- length(levels(df_16s_phylum$taxon))
pal_base <- c(RColorBrewer::brewer.pal(min(n_taxa, 12), "Paired"))
if (n_taxa > 12) pal_base <- c(pal_base, rep("#888888", n_taxa - 12))
names(pal_base) <- levels(df_16s_phylum$taxon)
pal_base["Other"] <- "#CCCCCC"
if ("Unclassified" %in% names(pal_base)) pal_base["Unclassified"] <- "#999999"

ggplot(df_16s_phylum, aes(x = site, y = mean_abund, fill = taxon)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = pal_base) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = "Phylum",
       title = "16S rRNA - Phylum-level composition") +
  theme_Publication() +
  guides(fill = guide_legend(ncol = 2))
ggsave("results/composition/16s_phylum_composition.pdf", width = 6, height = 6)

# Genus-level plot
df_16s_genus <- prep_abundance(df_16s, "Genus", top_n = 15, group_col = "site")

n_taxa_g <- length(levels(df_16s_genus$taxon))
pal_genus <- c(RColorBrewer::brewer.pal(min(n_taxa_g, 12), "Paired"))
if (n_taxa_g > 12) pal_genus <- c(pal_genus, colorRampPalette(c("#66C2A5", "#FC8D62", "#8DA0CB"))(n_taxa_g - 12))
names(pal_genus) <- levels(df_16s_genus$taxon)
pal_genus["Other"] <- "#CCCCCC"
if ("Unclassified" %in% names(pal_genus)) pal_genus["Unclassified"] <- "#999999"

ggplot(df_16s_genus, aes(x = site, y = mean_abund, fill = taxon)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = pal_genus) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = "Genus",
       title = "16S rRNA - Genus-level composition") +
  theme_Publication() +
  guides(fill = guide_legend(ncol = 2))
ggsave("results/composition/16s_genus_composition.pdf", width = 7, height = 7)

# ---- Shotgun composition plots ----
# Note: MetaPhlAn reports relative abundances at the species level only. Reads that
# cannot be classified to species (unclassified reads, higher-rank-only assignments,
# viral/eukaryotic reads) are excluded from species-level profiles. As a result,
# per-sample species abundances typically sum to 60-90% rather than 100%. The bars
# therefore reflect the classified bacterial/archaeal fraction, not the full community.
process_shotgun <- function(shotgun_obj, site_label) {
  counts_long <- shotgun_obj$counts |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "Species", values_to = "Abundance")

  # MetaPhlAn values are 0-100 scale, convert to proportions; NA = absent in one batch
  counts_long$Abundance <- replace_na(counts_long$Abundance, 0) / 100

  # Join taxonomy
  tax <- shotgun_obj$tax_table |> select(Species, Phylum, Genus)
  counts_long <- counts_long |> left_join(tax, by = "Species")

  # Join ethnicity
  sample_meta <- shotgun_obj$sample_data |>
    rownames_to_column("sample_id") |>
    select(sample_id, EthnicityTotal)
  counts_long <- counts_long |> left_join(sample_meta, by = "sample_id")

  # Remove samples without metadata
  counts_long <- counts_long |> filter(!is.na(EthnicityTotal))
  counts_long$site <- site_label
  counts_long
}

shotgun_throat <- readRDS("data/processed/shotgun_throat.RDS")
shotgun_tongue <- readRDS("data/processed/shotgun_tongue.RDS")

df_shotgun_throat <- process_shotgun(shotgun_throat, "Throat")
df_shotgun_tongue <- process_shotgun(shotgun_tongue, "Tongue")

# Shotgun throat - Phylum
df_sg_throat_phy <- prep_abundance(df_shotgun_throat, "Phylum", top_n = 10, group_col = "EthnicityTotal")

n_taxa_sg <- length(levels(df_sg_throat_phy$taxon))
pal_sg <- c(RColorBrewer::brewer.pal(min(n_taxa_sg, 12), "Paired"))
if (n_taxa_sg > 12) pal_sg <- c(pal_sg, rep("#888888", n_taxa_sg - 12))
names(pal_sg) <- levels(df_sg_throat_phy$taxon)
pal_sg["Other"] <- "#CCCCCC"
if ("Unclassified" %in% names(pal_sg)) pal_sg["Unclassified"] <- "#999999"

ggplot(df_sg_throat_phy, aes(x = EthnicityTotal, y = mean_abund, fill = taxon)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = pal_sg) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = "Phylum",
       title = "Shotgun metagenomics - Throat (Phylum)") +
  theme_Publication() +
  guides(fill = guide_legend(ncol = 2))
ggsave("results/composition/shotgun_phylum_throat.pdf", width = 6, height = 6)

# Shotgun tongue - Phylum
df_sg_tongue_phy <- prep_abundance(df_shotgun_tongue, "Phylum", top_n = 10, group_col = "EthnicityTotal")

pal_sg_tongue <- c(RColorBrewer::brewer.pal(min(length(levels(df_sg_tongue_phy$taxon)), 12), "Paired"))
if (length(levels(df_sg_tongue_phy$taxon)) > 12) pal_sg_tongue <- c(pal_sg_tongue, rep("#888888", length(levels(df_sg_tongue_phy$taxon)) - 12))
names(pal_sg_tongue) <- levels(df_sg_tongue_phy$taxon)
pal_sg_tongue["Other"] <- "#CCCCCC"
if ("Unclassified" %in% names(pal_sg_tongue)) pal_sg_tongue["Unclassified"] <- "#999999"

ggplot(df_sg_tongue_phy, aes(x = EthnicityTotal, y = mean_abund, fill = taxon)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = pal_sg_tongue) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = "Phylum",
       title = "Shotgun metagenomics - Tongue (Phylum)") +
  theme_Publication() +
  guides(fill = guide_legend(ncol = 2))
ggsave("results/composition/shotgun_phylum_tongue.pdf", width = 6, height = 6)

# Shotgun throat - Genus
df_sg_throat_gen <- prep_abundance(df_shotgun_throat, "Genus", top_n = 15, group_col = "EthnicityTotal")

n_taxa_sg_g <- length(levels(df_sg_throat_gen$taxon))
pal_sg_g <- c(RColorBrewer::brewer.pal(min(n_taxa_sg_g, 12), "Paired"))
if (n_taxa_sg_g > 12) pal_sg_g <- c(pal_sg_g, colorRampPalette(c("#66C2A5", "#FC8D62", "#8DA0CB"))(n_taxa_sg_g - 12))
names(pal_sg_g) <- levels(df_sg_throat_gen$taxon)
pal_sg_g["Other"] <- "#CCCCCC"
if ("Unclassified" %in% names(pal_sg_g)) pal_sg_g["Unclassified"] <- "#999999"

ggplot(df_sg_throat_gen, aes(x = EthnicityTotal, y = mean_abund, fill = taxon)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = pal_sg_g) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = "Genus",
       title = "Shotgun metagenomics - Throat (Genus)") +
  theme_Publication() +
  guides(fill = guide_legend(ncol = 3))
ggsave("results/composition/shotgun_genus_throat.pdf", width = 8, height = 7)

# Shotgun tongue - Genus
df_sg_tongue_gen <- prep_abundance(df_shotgun_tongue, "Genus", top_n = 15, group_col = "EthnicityTotal")

n_taxa_sg_tg <- length(levels(df_sg_tongue_gen$taxon))
pal_sg_tg <- c(RColorBrewer::brewer.pal(min(n_taxa_sg_tg, 12), "Paired"))
if (n_taxa_sg_tg > 12) pal_sg_tg <- c(pal_sg_tg, colorRampPalette(c("#66C2A5", "#FC8D62", "#8DA0CB"))(n_taxa_sg_tg - 12))
names(pal_sg_tg) <- levels(df_sg_tongue_gen$taxon)
pal_sg_tg["Other"] <- "#CCCCCC"
if ("Unclassified" %in% names(pal_sg_tg)) pal_sg_tg["Unclassified"] <- "#999999"

ggplot(df_sg_tongue_gen, aes(x = EthnicityTotal, y = mean_abund, fill = taxon)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = pal_sg_tg) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = "Genus",
       title = "Shotgun metagenomics - Tongue (Genus)") +
  theme_Publication() +
  guides(fill = guide_legend(ncol = 3))
ggsave("results/composition/shotgun_genus_tongue.pdf", width = 8, height = 7)
