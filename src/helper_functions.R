attribute_stats <- function(df, na.rm = TRUE) {
  
  ## Iterate over column names -------------------------------------------------
  col_list <- lapply(names(df), function(col) {
    vec <- df[[col]]
    
    ## basic descriptors -------------------------------------------------------
    is_num   <- is.numeric(vec)                   # numeric ? → mean/SD/hist
    cnt      <- length(vec)                       # total # of rows
    nulls    <- sum(is.na(vec))                    # NA count
    distinct_cnt <- length(unique(vec[!is.na(vec)]))  # distinct (ignores NAs)
    
    ## numeric‑only statistics --------------------------------------------------
    mn  <- if (is_num) min(vec, na.rm = TRUE) else NA_real_
    mx  <- if (is_num) max(vec, na.rm = TRUE) else NA_real_
    av  <- if (is_num) mean(vec, na.rm = na.rm) else NA_real_
    me <- if (is_num) median(vec, na.rm = na.rm) else NA_real_
    sdv <- if (is_num) sd(vec,   na.rm = na.rm) else NA_real_
    sov <- if (is_num) sum(vec,   na.rm = na.rm) else NA_real_ # sum of values
    sosv <- if (is_num) sum((as.numeric(vec))^2, na.rm = na.rm) else NA_real_ # sum of squred values
    n <- if (is_num) length(vec) - sum(is.na(vec)) else NA_real_# number of non NA values
    
    ## Return a **named list** that represents one row of the final table ------
    list(
      name        = col,
      count       = cnt,
      nulls       = nulls,
      distinct    = distinct_cnt,
      minimum     = mn,
      maximum     = mx,
      average     = av,
      median      = me,
      stddev      = sdv,
      sumofvals   = sov,
      sumofsqvals = sosv,
      n           = n
    )
  })}

las_info <- function(las){
  info <- list()
  info$bbox <- as.list(lidR::st_bbox(las))
  info$point_count <- nrow(las@data)
  info$crs <- jsonlite::toJSON(lapply(lidR::st_crs(las), as.list))
  info$dimensions <- attribute_stats(las@data)
  return(info)
}

multi_tile_variance <- function(n, sov, sosv) {
  # n    : numeric vector of non‑NA counts per tile
  # sov  : numeric vector of sum of values per tile
  # sosv : numeric vector of sum of squared values per tile
  
  N <- sum(n)
  if (N == 0) {
    return(list(mean = NA_real_, variance = NA_real_, stddev = NA_real_))
  }
  
  mv  <- sum(sov) / N
  var <- (sum(sosv) / N) - mv^2
  stdd <- sqrt(var)
  
  list(mean = mv, variance = var, stddev = stdd)
}

chull_to_polygon <- function(coords, crs = sf::NA_crs_) {
  stopifnot(ncol(coords) == 2)
  
  m <- as.matrix(coords)
  # close ring if needed
  if (!all(m[1, ] == m[nrow(m), ])) {
    m <- rbind(m, m[1, ])
  }
  
  sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(m)), crs = crs)
  )
}

polygons_to_multipolygon <- function(polys, files = NULL) {
  # polys: list of sf objects with POLYGON geometry
  geoms <- do.call(c, lapply(polys, sf::st_geometry))  # sfc_POLYGON
  crs   <- sf::st_crs(geoms)
  
  # single MULTIPOLYGON geometry combining all polygons
  mp <- sf::st_cast(sf::st_combine(geoms), "MULTIPOLYGON")
  mp <- sf::st_sfc(mp, crs = crs)
  
  if (is.null(files)) {
    return(sf::st_sf(geometry = mp))
  }
  
  # attach filenames as a list-column (one vector per MULTIPOLYGON feature)
  if (length(files) != length(geoms)) {
    stop("Length of 'files' must match number of input polygons.")
  }
  
  sf::st_sf(
    files    = I(list(files)),
    geometry = mp
  )
}

# Is the CRS defined / usable by sf?
is_crs_defined <- function(crs) {
  !is.na(crs) && !is.na(crs$wkt)
}

# Can this CRS be transformed to some target (e.g. EPSG:4326)?
is_crs_transformable <- function(crs, target = 4326) {
  if (!is_crs_defined(crs)) return(FALSE)
  
  # create a dummy geometry with this CRS and try to transform it
  g <- sf::st_sfc(sf::st_point(c(0, 0)), crs = crs)
  res <- try(sf::st_transform(g, target), silent = TRUE)
  
  !inherits(res, "try-error")
}

is_valid_wgs84_point <- function(lon, lat) {
  lon <- suppressWarnings(as.numeric(lon))[1]
  lat <- suppressWarnings(as.numeric(lat))[1]

  is.finite(lon) && is.finite(lat) &&
    lon >= -180 && lon <= 180 &&
    lat >= -90 && lat <= 90
}

unavailable_centroid <- function(source = "unavailable") {
  list(
    longitude = NA_real_,
    latitude = NA_real_,
    crs = "EPSG:4326",
    source = source,
    available = FALSE
  )
}

collection_centroid <- function(mp_sf, geometry_is_wgs84, approx_lon = NA_real_, approx_lat = NA_real_) {
  if (is_valid_wgs84_point(approx_lon, approx_lat)) {
    return(list(
      longitude = as.numeric(approx_lon)[1],
      latitude = as.numeric(approx_lat)[1],
      crs = "EPSG:4326",
      source = "approximate_input",
      available = TRUE
    ))
  }

  if (isTRUE(geometry_is_wgs84)) {
    centroid <- tryCatch({
      point <- suppressWarnings(sf::st_centroid(mp_sf))
      coords <- sf::st_coordinates(point)[1, ]

      if (is_valid_wgs84_point(coords[["X"]], coords[["Y"]])) {
        list(
          longitude = unname(coords[["X"]]),
          latitude = unname(coords[["Y"]]),
          crs = "EPSG:4326",
          source = "geometry",
          available = TRUE
        )
      } else {
        NULL
      }
    }, error = function(e) NULL)

    if (is.null(centroid)) {
      centroid <- tryCatch({
        bbox <- sf::st_bbox(mp_sf)
        lon <- mean(c(bbox[["xmin"]], bbox[["xmax"]]))
        lat <- mean(c(bbox[["ymin"]], bbox[["ymax"]]))

        if (is_valid_wgs84_point(lon, lat)) {
          list(
            longitude = unname(lon),
            latitude = unname(lat),
            crs = "EPSG:4326",
            source = "geometry_bbox",
            available = TRUE
          )
        } else {
          NULL
        }
      }, error = function(e) NULL)
    }

    if (!is.null(centroid)) {
      return(centroid)
    }
  }

  unavailable_centroid("unavailable")
}
