source("src/helper_functions.R")

centroid <- collection_centroid(
  mp_sf = stop("geometry should not be evaluated when approximate coordinates are valid"),
  geometry_is_wgs84 = TRUE,
  approx_lon = 7.85,
  approx_lat = 47.99
)

expected <- list(
  longitude = 7.85,
  latitude = 47.99,
  crs = "EPSG:4326",
  source = "approximate_input",
  available = TRUE
)

if (!identical(centroid, expected)) {
  stop(
    "Expected approximate coordinates to override geometry centroid, got: ",
    jsonlite::toJSON(centroid, auto_unbox = TRUE, na = "null")
  )
}
