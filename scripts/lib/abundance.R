## Taxon labels and relative-abundance inputs shared by DA reports.
## Callers choose significance thresholds and features before preparing plots.
make_taxon_labels <- function(results) {
    results |>
        dplyr::mutate(label = dplyr::case_when(
            !is.na(Tax) & Tax != "" ~ Tax,
            TRUE ~ feature
        )) |>
        dplyr::group_by(label) |>
        dplyr::mutate(label = if (dplyr::n() > 1) paste0(label, " (", feature, ")") else label) |>
        dplyr::ungroup()
}

mean_abundance_by_ethnicity <- function(ps, meta, features) {
    ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    abund_mat <- as.data.frame(phyloseq::otu_table(ps_rel))
    if (phyloseq::taxa_are_rows(ps_rel)) {
        abund_mat <- abund_mat[features, , drop = FALSE]
    } else {
        abund_mat <- t(abund_mat)[features, , drop = FALSE]
    }

    ## Compute mean abundance per ethnicity. Both sapply and vapply
    ## silently simplify a length-1 per-group result down to a bare
    ## vector (losing the row dimension) when only one taxon is
    ## significant, so build the matrix explicitly via cbind instead,
    ## which always preserves it regardless of row count.
    mean_abund <- sapply(levels(meta$EthnicityTotal), function(eth) {
        samples <- rownames(meta[meta$EthnicityTotal == eth, , drop = FALSE])
        rowMeans(abund_mat[, samples, drop = FALSE])
    }, simplify = FALSE) |> do.call(cbind, args = _)

    mean_abund
}

prepare_taxon_boxplot_data <- function(ps, meta, top_box) {
    ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    abund_df <- as.data.frame(phyloseq::otu_table(ps_rel))
    if (phyloseq::taxa_are_rows(ps_rel)) {
        abund_df <- as.data.frame(t(abund_df))
    }

    ## Select top taxa and pivot
    box_df <- abund_df[, top_box$feature, drop = FALSE] |>
        tibble::rownames_to_column("sample_id") |>
        tidyr::pivot_longer(-sample_id, names_to = "feature", values_to = "rel_abund") |>
        dplyr::left_join(
            meta |> tibble::rownames_to_column("sample_id") |>
                dplyr::select(sample_id, EthnicityTotal),
            by = "sample_id"
        ) |>
        dplyr::left_join(top_box |> dplyr::select(feature, label, qval), by = "feature") |>
        dplyr::mutate(label = paste0(label, "\n(q = ",
                              formatC(qval, format = "e", digits = 1), ")"))

    box_df
}
