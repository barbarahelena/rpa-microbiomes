## Shared publication theme and ethnicity palette.

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

ethnicity_colours <- c(
    "Dutch"                  = "#1F78B4",
    "South-Asian Surinamese" = "#E31A1C",
    "African Surinamese"     = "#33A02C",
    "Javanese Surinamese"    = "#6A3D9A",
    "Other"                  = "#B15928",
    "Ghanaian"               = "#FF7F00",
    "Turkish"                = "#E7298A",
    "Moroccan"               = "#D4AC0D"
)

## Draw an already-computed p-value at the original ggpubr panel position.
## Positioning happens after scale training, so annotations do not change
## the data range used to determine significance-bracket tip lengths.
stat_cached_pvalue <- function(label, label.y) {
    ggplot2::layer(
        stat = ggplot2::ggproto(NULL, ggplot2::Stat,
            required_aes = c("x", "y"),
            compute_panel = function(data, scales, label, label.y) {
                data.frame(x = scales$x$dimension()[1], y = label.y,
                           label = label, hjust = 0.2, vjust = 0.5)
            }),
        geom = "text", position = "identity", inherit.aes = TRUE,
        params = list(label = label, label.y = label.y, vjust = 0)
    )
}
