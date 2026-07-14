## Overlap of significant differential-abundance taxa across ethnicity pairs
## (16S throat and nose)
##
## Reads the per-pair MaAsLin2 results already written by
## 8_differential_abundance_16s.R (maaslin2_ethnicity_results_16s_<site>_
## <pair>.csv) and builds an UpSet plot per site showing which ASVs are
## significant (q < 0.05) in which pairwise ethnicity comparisons, plus a CSV
## of the ASVs shared by 2+ pairs for follow-up. Doesn't refit anything, so
## it's cheap to re-run while iterating on the plot.

## Libraries
library(here)
library(tidyverse)
library(UpSetR)

## Setup
setwd(here::here())

## Use the same outdir as script 8's DIFF_AB_TEST_N, so this can be pointed
## at a test-mode run's output for a quick check.
test_n <- suppressWarnings(as.integer(Sys.getenv("DIFF_AB_TEST_N", "")))
maaslin_dir <- if (!is.na(test_n)) "results/differential_abundance_test/maaslin2" else "results/differential_abundance/maaslin2"
outdir <- if (!is.na(test_n)) "results/differential_abundance_test/upset" else "results/differential_abundance/upset"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Significance threshold for set membership - matches the "strict" q < 0.05
## cutoff script 8 already uses for its forest plots and top-taxa boxplots.
SIG_THRESHOLD <- 0.05

sites <- c("throat", "nose")

for (site_name in sites) {
    files <- list.files(
        maaslin_dir,
        pattern = paste0("^maaslin2_ethnicity_results_16s_", site_name, "_.+\\.csv$"),
        full.names = TRUE
    )

    if (length(files) == 0) {
        cat(site_name, "- no MaAsLin2 pair results found in", maaslin_dir,
            "- run 8_differential_abundance_16s.R first. Skipping.\n\n")
        next
    }

    ## pair_tag (e.g. "Dutch_vs_SAS") is everything after the site name
    pair_tags <- files |>
        basename() |>
        str_remove(paste0("^maaslin2_ethnicity_results_16s_", site_name, "_")) |>
        str_remove("\\.csv$")

    results <- map(files, read_csv, show_col_types = FALSE) |>
        set_names(pair_tags)

    ## Taxonomy lookup (feature -> readable label), same for every pair at a
    ## site since they share the same ASV table
    tax_lookup <- bind_rows(results) |>
        distinct(feature, Tax)

    ## Named list of significant ASVs per pair, for UpSetR::fromList()
    sig_sets <- map(results, ~ .x |> filter(qval < SIG_THRESHOLD) |> pull(feature))
    sig_sets <- sig_sets[lengths(sig_sets) > 0]

    if (length(sig_sets) < 2) {
        cat(site_name, "- fewer than 2 pairs have any significant taxa",
            "(q <", SIG_THRESHOLD, ") - skipping UpSet plot.\n\n")
        next
    }

    ## ---- UpSet plot ----
    membership <- UpSetR::fromList(sig_sets)
    pdf(file.path(outdir, paste0("upset_diffabund_16s_", site_name, ".pdf")),
        width = 9, height = 6, onefile = FALSE)
    print(UpSetR::upset(
        membership,
        sets = names(sig_sets),
        nintersects = 30,
        order.by = "freq",
        mainbar.y.label = paste0("Shared significant ASVs (q < ", SIG_THRESHOLD, ")"),
        sets.x.label = "Significant ASVs per pair",
        text.scale = 1.1
    ))
    dev.off()

    ## ---- Overlap detail: ASVs significant in 2+ pairs, with taxonomy ----
    ## fromList() doesn't keep the original element names as rownames - it
    ## just returns 1:N in the order of unique(unlist(sig_sets)), so rebuild
    ## that same vector to recover which feature each row corresponds to.
    membership_df <- membership |>
        mutate(feature = unique(unlist(sig_sets))) |>
        rowwise() |>
        mutate(n_pairs = sum(c_across(all_of(names(sig_sets)))),
               pairs = paste(names(sig_sets)[c_across(all_of(names(sig_sets))) == 1],
                            collapse = "; ")) |>
        ungroup() |>
        filter(n_pairs >= 2) |>
        left_join(tax_lookup, by = "feature") |>
        select(feature, Tax, n_pairs, pairs) |>
        arrange(desc(n_pairs))

    write_csv(membership_df,
              file.path(outdir, paste0("upset_overlap_16s_", site_name, ".csv")))

    cat("Completed:", site_name, "-", length(sig_sets), "pairs with signal,",
        nrow(membership_df), "ASVs shared by 2+ pairs\n\n")
}
