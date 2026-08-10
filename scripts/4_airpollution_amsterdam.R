
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
  st_simplify(dTolerance = 2) |>
  st_make_valid() |>
  st_transform(4326)

amsterdam_boundary_wgs84 <- st_transform(amsterdam_boundary, 4326)

pollutant_cols <- c(
  "PM10 2013" = "conc_ALO_pm10_2013", "PM10 2014" = "conc_ALO_pm10_2014", "PM10 2015" = "conc_ALO_pm10_2015",
  "PM2.5 2013" = "conc_ALO_pm25_2013", "PM2.5 2014" = "conc_ALO_pm25_2014", "PM2.5 2015" = "conc_ALO_pm25_2015",
  "NO2 2014" = "conc_ALO_no2_2014", "NO2 2015" = "conc_ALO_no2_2015",
  "EC 2013" = "conc_ALO_ec_2013", "EC 2014" = "conc_ALO_ec_2014", "EC 2015" = "conc_ALO_ec_2015"
)
default_pollutant <- "PM2.5 2015"

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
