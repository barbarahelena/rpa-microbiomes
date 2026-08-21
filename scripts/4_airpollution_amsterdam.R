
## Data cleaning of air pollution data

## Libraries
library(here)
library(tidyverse)
library(sf)
library(leaflet)
library(htmlwidgets)
library(htmltools)

setwd(here::here())
dir.create("results/airpollution", recursive = TRUE, showWarnings = FALSE)

# Data
meta1 <- haven::read_sav("data/raw/250606_HELIUS data Barbara Verhaar.sav")
meta2 <- haven::read_sav("data/raw/231108_HELIUS data Barbara Verhaar_GECCO.sav")
meta2 <- meta2 |> dplyr::select(ID, conc_ALO_pm10_2013:conc_ALO_ec_2015)

pc6 <- read.csv("data/raw/PC6_2022_ALO_2013_2015.csv")

## Map of Amsterdam air pollution by PC6 postcode
# PC6 geometries: CBS Postcode6 2022 (PDOK), downloaded to data/raw/geo/cbs_pc6_2022.gpkg
# Amsterdam boundary: CBS Gebiedsindelingen 2022, gemeente_gegeneraliseerd layer,
#   pre-filtered to Amsterdam and saved to data/raw/geo/amsterdam_gemeente.gpkg

amsterdam_boundary <- st_read("data/raw/geo/amsterdam_gemeente.gpkg", quiet = TRUE)

# Spatial pre-filter (bbox, uses gpkg spatial index) then exact intersect with the boundary
bbox_wkt <- st_as_text(st_as_sfc(st_bbox(amsterdam_boundary)))

pc6_geo <- st_read(
  "data/raw/geo/cbs_pc6_2022.gpkg",
  wkt_filter = bbox_wkt,
  quiet = TRUE
)

pc6_geo <- pc6_geo[st_intersects(pc6_geo, amsterdam_boundary, sparse = FALSE)[, 1], ]

pc6_amsterdam <- pc6_geo |>
  select(postcode6, geom) |>
  left_join(pc6, by = "postcode6") |>
  ## RIVM/ALO background concentrations are never literally 0 in an urban
  ## Dutch setting - a stored 0 is a "no estimate for this PC6" placeholder,
  ## not a real reading (e.g. postcode6 6245ES/4341RM/etc. are 0 across all
  ## three years for every pollutant). Left as 0 these drag both the
  ## multi-year averages and the map's colour scale down. Recode to NA so
  ## rowMeans(na.rm = TRUE) below skips them like any other missing value.
  mutate(across(starts_with("conc_ALO_"), ~ na_if(.x, 0))) |>
  st_simplify(dTolerance = 2) |>
  st_make_valid() |>
  st_transform(4326)

## Drop a handful of corrupted source polygons: self-intersecting features
## whose own bounding-box diagonal (in degrees, post-transform) is wildly
## larger than a real PC6 postcode (postcodes are at most a few hundred
## metres across). Left in, a single such polygon (e.g. postcode 3053BM,
## actually in Rotterdam, ~0.57 degrees across vs ~0.001-0.04 for every
## legitimate postcode) blows out the fixed-aspect map's bounding box and
## squeezes Amsterdam into a corner with the rest left blank. Threshold
## (0.05 deg, ~5km) is well above the legitimate 99.9th percentile
## (~0.038 deg) and well below the corrupted outlier.
bbox_diag <- function(x) {
  b <- st_bbox(x)
  sqrt((b["xmax"] - b["xmin"])^2 + (b["ymax"] - b["ymin"])^2)
}
feature_diag <- vapply(st_geometry(pc6_amsterdam), bbox_diag, numeric(1))
is_corrupted <- !is.na(feature_diag) & feature_diag > 0.05
if (any(is_corrupted)) {
  cat("Dropping", sum(is_corrupted), "corrupted PC6 polygon(s) with implausible extent:",
      paste(pc6_amsterdam$postcode6[is_corrupted], collapse = ", "), "\n")
}
pc6_amsterdam <- pc6_amsterdam[!is_corrupted, ]

pc6_amsterdam <- pc6_amsterdam |>
  mutate(
    pm10_avg_2013_2015 = rowMeans(
      across(c(conc_ALO_pm10_2013, conc_ALO_pm10_2014, conc_ALO_pm10_2015)),
      na.rm = TRUE
    ),
    pm25_avg_2013_2015 = rowMeans(
      across(c(conc_ALO_pm25_2013, conc_ALO_pm25_2014, conc_ALO_pm25_2015)),
      na.rm = TRUE
    ),
    ## NO2 is only available for 2014-2015 (see pollutant_cols below)
    no2_avg_2014_2015 = rowMeans(
      across(c(conc_ALO_no2_2014, conc_ALO_no2_2015)),
      na.rm = TRUE
    ),
    ec_avg_2013_2015 = rowMeans(
      across(c(conc_ALO_ec_2013, conc_ALO_ec_2014, conc_ALO_ec_2015)),
      na.rm = TRUE
    )
  )

amsterdam_boundary_wgs84 <- st_transform(amsterdam_boundary, 4326)

## Cache the processed geometries (PC6 polygons + Amsterdam boundary, WGS84)
## so downstream scripts (Figure 1) can build a small static map panel
## without repeating this spatial read/join.
saveRDS(list(pc6 = pc6_amsterdam, boundary = amsterdam_boundary_wgs84),
        "results/airpollution/amsterdam_pc6_geo.rds")

## Static, non-interactive maps (one per pollutant, multi-year mean) -
## standalone sanity-check exports; the small Figure 1 panel is rebuilt from
## the cached geometries above (see 12_figure1_assembly.R) rather than
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
