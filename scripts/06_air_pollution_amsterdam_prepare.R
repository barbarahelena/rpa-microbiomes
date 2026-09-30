## Prepare Amsterdam postcode geometries and air-pollution estimates.
## 1. Select postcode polygons intersecting Amsterdam and join exposure data.
## 2. Treat zero estimates as missing, simplify/repair geometries, and remove
##    polygons with implausible extents.
## 3. Calculate multi-year pollutant means and cache WGS84 map inputs.

## Libraries
library(here)
library(tidyverse)
library(sf)
library(leaflet)
library(htmlwidgets)
library(htmltools)

setwd(here::here())
dir.create("results/airpollution", recursive = TRUE, showWarnings = FALSE)

## ---- Load exposures and select Amsterdam postcode polygons ----
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

## ---- Join exposures and clean estimates and geometries ----
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

## ---- Calculate multi-year pollutant means ----
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

## ---- Cache map inputs for reports and Figure 1 ----
## Cache the processed geometries (PC6 polygons + Amsterdam boundary, WGS84)
## so downstream scripts (Figure 1) can build a small static map panel
## without repeating this spatial read/join.
saveRDS(list(pc6 = pc6_amsterdam, boundary = amsterdam_boundary_wgs84),
        "results/airpollution/amsterdam_pc6_geo.rds")

