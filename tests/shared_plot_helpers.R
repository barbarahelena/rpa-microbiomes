## Small synthetic checks for reusable plotting inputs; no cohort data or models.
library(tidyverse)
library(phyloseq)
source("scripts/lib/plot_annotations.R")
source("scripts/lib/abundance.R")
source("scripts/lib/diversity_plots.R")

## Brackets retain per-facet ranges, strict thresholds, ordering and caps.
pairs <- tibble(metric = rep(c("a", "b"), each = 3),
                group1 = "Dutch", group2 = "Turkish",
                p.adj = c(0.04, 0.001, 0.05, 0.00001, NA, 0.01),
                min_val = 0, max_val = rep(c(10, 100), each = 3))
brackets <- pairs |> group_by(metric) |>
    prepare_significance_brackets(spacing = 0.1, max_pairs = 1) |> ungroup() |> arrange(p.adj)
stopifnot(identical(brackets$metric, c("b", "a")),
          identical(brackets$p.adj.label, c("****", "**")),
          identical(brackets$y.position, c(110, 11)))
empty <- pairs |> filter(is.na(p.adj) | p.adj >= 0.05) |>
    prepare_significance_brackets(spacing = 0.06)
stopifnot(nrow(empty) == 0L, is.numeric(empty$y.position))
pm <- matrix(c(NA, 0.02, NA, NA), 2, dimnames = list(c("a", "b"), c("a", "b")))
tidy <- tidy_pairwise_pvalues(pm)
stopifnot(nrow(tidy) == 1L, as.character(tidy$group1) == "b",
          as.character(tidy$group2) == "a", tidy$p.adj == 0.02)

## Labels fall back to IDs and distinguish repeated names without reordering.
taxa <- tibble(feature = c("f1", "f2", "f3", "f4"),
               Tax = c("Shared", "Shared", "", NA_character_), qval = 0.01)
labelled <- make_taxon_labels(taxa)
stopifnot(identical(labelled$feature, taxa$feature),
          identical(labelled$label, c("Shared (f1)", "Shared (f2)", "f3", "f4")))

## Abundance helpers preserve one-taxon dimensions and align by sample name,
## including when metadata order differs or taxa occupy columns instead of rows.
counts <- matrix(c(1, 3, 2, 2, 3, 1, 4, 0), nrow = 2,
                 dimnames = list(c("f1", "f2"), paste0("s", 1:4)))
meta <- data.frame(EthnicityTotal = factor(c("Dutch", "Dutch", "Turkish", "Turkish")),
                   row.names = colnames(counts))
ps <- phyloseq(otu_table(counts, taxa_are_rows = TRUE), sample_data(meta))
ps_transposed <- phyloseq(otu_table(t(counts), taxa_are_rows = FALSE), sample_data(meta))
shuffled <- meta[c(4, 1, 3, 2), , drop = FALSE]
a <- mean_abundance_by_ethnicity(ps, shuffled, "f1")
b <- mean_abundance_by_ethnicity(ps_transposed, shuffled, "f1")
stopifnot(identical(dim(a), c(1L, 2L)), identical(a, b),
          identical(as.numeric(a), c(0.375, 0.875)))
top <- labelled[1, ]
a <- prepare_taxon_boxplot_data(ps, shuffled, top)
b <- prepare_taxon_boxplot_data(ps_transposed, shuffled, top)
stopifnot(identical(a, b), identical(a$sample_id, colnames(counts)),
          identical(a$rel_abund, c(0.25, 0.5, 0.75, 1)),
          identical(a$EthnicityTotal, meta$EthnicityTotal))

## Centroids are added only for four or more observed groups.
set.seed(7)
for (n_groups in c(3, 4)) {
    metadata <- data.frame(EthnicityTotal = factor(rep(LETTERS[seq_len(n_groups)], each = 10)))
    pcoa <- list(values = list(Eigenvalues = c(6, 3, 1)),
                 vectors = matrix(rnorm(nrow(metadata) * 2), ncol = 2))
    prepared <- prepare_pcoa_plot_data(metadata, pcoa)
    palette <- setNames(c("red", "blue", "green", "orange")[seq_len(n_groups)], levels(metadata$EthnicityTotal))
    plot <- plot_ethnicity_pcoa(prepared$ord_df, prepared$n_groups, palette)
    built <- ggplot_build(plot)
    stopifnot(identical(prepared$var_explained, c(60, 30, 10)),
              length(built$data) == if (n_groups > 3) 3L else 2L)
    if (n_groups > 3) stopifnot(nrow(built$data[[3]]) == n_groups)
}
cat("Shared plot helper edge cases passed.\n")
