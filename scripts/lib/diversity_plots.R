## Base plots shared by standalone reports and manuscript panels.
## These functions consume saved analysis data and never fit models.
prepare_pcoa_plot_data <- function(meta, pcoa) {
    eig <- pcoa$values$Eigenvalues
    list(
        var_explained = round(100 * eig / sum(eig), 1),
        ord_df = data.frame(PCo1 = pcoa$vectors[, 1],
                            PCo2 = pcoa$vectors[, 2],
                            EthnicityTotal = meta$EthnicityTotal),
        n_groups = nlevels(droplevels(meta$EthnicityTotal))
    )
}

plot_ethnicity_pcoa <- function(ord_df, n_groups, eth_colours) {
    p <- ggplot2::ggplot(ord_df, ggplot2::aes(x = PCo1, y = PCo2, colour = EthnicityTotal)) +
        ggplot2::geom_point(alpha = 0.5, size = 1) +
        ggplot2::stat_ellipse(level = 0.95, linewidth = 0.8)
    if (n_groups > 3) {
        centroids <- ord_df |>
            dplyr::group_by(EthnicityTotal) |>
            dplyr::summarise(PCo1 = mean(PCo1), PCo2 = mean(PCo2), .groups = "drop")
        p <- p + ggplot2::geom_point(data = centroids,
            ggplot2::aes(x = PCo1, y = PCo2, fill = EthnicityTotal),
            shape = 21, colour = "black", size = 4, stroke = 0.8) +
            ggplot2::scale_fill_manual(values = eth_colours, guide = "none")
    }
    p
}

plot_sampling_seasonality <- function(date_df, eth_colours, legend_title = ggplot2::waiver()) {
    month_starts <- lubridate::yday(as.Date(paste0("2001-", 1:12, "-01")))
    ggplot2::ggplot(date_df, ggplot2::aes(x = yday, colour = EthnicityTotal, fill = EthnicityTotal)) +
        ggplot2::geom_density(alpha = 0.15, linewidth = 0.8) +
        ggplot2::scale_colour_manual(values = eth_colours, name = legend_title) +
        ggplot2::scale_fill_manual(values = eth_colours, name = legend_title) +
        ggplot2::scale_x_continuous(breaks = month_starts, labels = month.abb,
                                  limits = c(1, 366), expand = c(0, 0))
}
