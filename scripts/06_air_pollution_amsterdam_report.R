## Report Amsterdam-wide air-pollution maps.
## 1. Load postcode exposures and boundaries from the preparation cache.
## 2. Export static PDF/PNG maps for PM10, PM2.5, NO2, and soot/EC.
## 3. Save an interactive HTML map with a pollutant selector.

library(here)
library(tidyverse)
library(sf)
library(leaflet)
library(htmlwidgets)
library(htmltools)

setwd(here::here())
dir.create("results/airpollution", recursive = TRUE, showWarnings = FALSE)

## ---- Load cached postcode exposures and boundaries ----
cache_path <- "results/airpollution/amsterdam_pc6_geo.rds"
if (!file.exists(cache_path)) stop("Missing cache: ", cache_path, "; run pixi run airpollution-amsterdam-prepare first.")
cache <- readRDS(cache_path)
pc6_amsterdam <- cache$pc6
amsterdam_boundary_wgs84 <- cache$boundary

## ---- Export static maps for each pollutant ----
## Static, non-interactive maps (one per pollutant, multi-year mean) -
## standalone sanity-check exports; the small Figure 1 panel is rebuilt from
## the cached geometries above (see 14_figure1.R) rather than
## sourcing these plots directly. coord_sf(expand = FALSE) plus zero plot
## margins keep the map filling its panel instead of floating in whitespace.
static_map_specs <- list(
  list(file = "pm10_2013_2015", col = "pm10_avg_2013_2015", title = "PM10 (2013-2015 mean)", legend = "PM10\n(µg/m³)"),
  list(file = "pm25_2013_2015", col = "pm25_avg_2013_2015", title = "PM2.5 (2013-2015 mean)", legend = "PM2.5\n(µg/m³)"),
  list(file = "no2_2014_2015",  col = "no2_avg_2014_2015",  title = "NO2 (2014-2015 mean)",   legend = "NO2\n(µg/m³)"),
  list(file = "ec_2013_2015",   col = "ec_avg_2013_2015",   title = "Soot/EC (2013-2015 mean)", legend = "EC\n(µg/m³)")
)

for (spec in static_map_specs) {
  static_map <- ggplot(pc6_amsterdam) +
    geom_sf(aes(fill = .data[[spec$col]]), colour = NA) +
    geom_sf(data = amsterdam_boundary_wgs84, fill = NA, colour = "black", linewidth = 0.4) +
    coord_sf(expand = FALSE) +
    scale_fill_viridis_c(name = spec$legend, option = "magma", direction = -1,
                          na.value = "grey85",
                          guide = guide_colorbar(barwidth = unit(0.3, "cm"), barheight = unit(3, "cm"))) +
    labs(title = paste0("Amsterdam - ", spec$title)) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.margin = margin(2, 2, 2, 2),
          legend.position = "right")

  ggsave(paste0("results/airpollution/amsterdam_", spec$file, "_static.pdf"), static_map, width = 6, height = 5)
  ggsave(paste0("results/airpollution/amsterdam_", spec$file, "_static.png"), static_map, width = 6, height = 5, dpi = 300)
}

## ---- Build and save the interactive pollutant map ----
pollutant_cols <- c(
  "PM10 2013" = "conc_ALO_pm10_2013", "PM10 2014" = "conc_ALO_pm10_2014", "PM10 2015" = "conc_ALO_pm10_2015",
  "PM10 avg 2013-2015" = "pm10_avg_2013_2015",
  "PM2.5 2013" = "conc_ALO_pm25_2013", "PM2.5 2014" = "conc_ALO_pm25_2014", "PM2.5 2015" = "conc_ALO_pm25_2015",
  "PM2.5 avg 2013-2015" = "pm25_avg_2013_2015",
  "NO2 2014" = "conc_ALO_no2_2014", "NO2 2015" = "conc_ALO_no2_2015",
  "NO2 avg 2014-2015" = "no2_avg_2014_2015",
  "EC 2013" = "conc_ALO_ec_2013", "EC 2014" = "conc_ALO_ec_2014", "EC 2015" = "conc_ALO_ec_2015",
  "EC avg 2013-2015" = "ec_avg_2013_2015"
)
default_pollutant <- "PM2.5 avg 2013-2015"

# Draw the postcode polygons once; pollutant switching recolors this single
# layer via JS instead of adding a duplicate geometry layer per pollutant
# (which previously produced an 11x larger, 130MB+ HTML file).
na_color <- "#cccccc"

colors_by_pollutant <- lapply(pollutant_cols, function(col) {
  values <- pc6_amsterdam[[col]]
  pal <- colorNumeric("YlOrRd", domain = values, na.color = na_color)
  setNames(as.list(pal(values)), pc6_amsterdam$postcode6)
}) |> setNames(names(pollutant_cols))

legend_html_by_pollutant <- lapply(names(pollutant_cols), function(label) {
  values <- pc6_amsterdam[[pollutant_cols[[label]]]]
  rng <- range(values, na.rm = TRUE)
  stops <- colorNumeric("YlOrRd", domain = values)(seq(rng[1], rng[2], length.out = 5))
  as.character(tags$div(
    style = "font: 12px sans-serif;",
    tags$b(label), tags$br(),
    tags$div(style = sprintf(
      "width:150px;height:12px;background:linear-gradient(to right, %s);",
      paste(stops, collapse = ", ")
    )),
    tags$div(
      style = "display:flex;justify-content:space-between;width:150px;",
      tags$span(round(rng[1], 1)), tags$span(round(rng[2], 1))
    )
  ))
}) |> setNames(names(pollutant_cols))

dropdown <- tags$div(
  style = "background:white;padding:6px;border-radius:4px;box-shadow:0 1px 4px rgba(0,0,0,0.3);",
  tags$select(
    id = "pollutant-select",
    lapply(names(pollutant_cols), function(nm) {
      tags$option(value = nm, nm, selected = if (nm == default_pollutant) "selected" else NULL)
    })
  ),
  tags$div(id = "pollutant-legend")
)

map <- leaflet(pc6_amsterdam) |>
  addProviderTiles("CartoDB.Positron") |>
  addPolygons(
    layerId = ~postcode6,
    fillColor = unlist(colors_by_pollutant[[default_pollutant]]),
    fillOpacity = 0.75,
    color = "white",
    weight = 0.3,
    label = ~postcode6
  ) |>
  addPolygons(data = amsterdam_boundary_wgs84, color = "black", weight = 2, fill = FALSE) |>
  addControl(html = dropdown, position = "topright")

map <- onRender(
  map,
  "
  function(el, x, data) {
    var map = this;
    document.getElementById('pollutant-legend').innerHTML = data.legends[data.default];
    document.getElementById('pollutant-select').addEventListener('change', function() {
      var pollutant = this.value;
      var colors = data.colors[pollutant];
      map.eachLayer(function(layer) {
        var id = layer.options && layer.options.layerId;
        if (id && colors[id]) { layer.setStyle({fillColor: colors[id]}); }
      });
      document.getElementById('pollutant-legend').innerHTML = data.legends[pollutant];
    });
  }
  ",
  data = list(colors = colors_by_pollutant, legends = legend_html_by_pollutant, default = default_pollutant)
)

saveWidget(map, "results/airpollution/amsterdam_pollution_map.html", selfcontained = TRUE)
