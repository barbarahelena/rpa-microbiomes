## Presentation helpers for already-computed pairwise tests.
## Input to prepare_significance_brackets has p.adj, group1/group2 and
## max_val/min_val columns. Existing grouping (e.g. metric) is preserved.
tidy_pairwise_pvalues <- function(p_values) {
    as.data.frame(as.table(p_values)) |>
        dplyr::filter(!is.na(Freq)) |>
        dplyr::rename(group1 = Var1, group2 = Var2, p.adj = Freq)
}

prepare_significance_brackets <- function(pairs, spacing, max_pairs = Inf) {
    pairs <- pairs |>
        dplyr::filter(p.adj < 0.05) |>
        dplyr::arrange(p.adj)
    if (is.finite(max_pairs)) pairs <- dplyr::slice_head(pairs, n = max_pairs)
    pairs |>
        dplyr::mutate(
            group1 = as.character(group1),
            group2 = as.character(group2),
            step = (max_val - min_val) * spacing,
            y.position = max_val + step * dplyr::row_number(),
            p.adj.label = dplyr::case_when(
                p.adj < 0.0001 ~ "****",
                p.adj < 0.001 ~ "***",
                p.adj < 0.01 ~ "**",
                TRUE ~ "*"
            )
        )
}
