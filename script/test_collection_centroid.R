source("src/helper_functions.R")
source("src/parameters.R")

parsed_location <- parse_approx_location("(7.85, 47.99)")
if (!identical(parsed_location, list(lon = 7.85, lat = 47.99))) {
  stop(
    "Expected approximate location tuple to parse as lon/lat, got: ",
    jsonlite::toJSON(parsed_location, auto_unbox = TRUE, na = "null")
  )
}

centroid <- collection_centroid(
  mp_sf = stop("geometry should not be evaluated when approximate coordinates are valid"),
  geometry_is_wgs84 = TRUE,
  approx_lon = parsed_location$lon,
  approx_lat = parsed_location$lat
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
