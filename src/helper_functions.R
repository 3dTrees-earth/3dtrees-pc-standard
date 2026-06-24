attribute_stats <- function(df, na.rm = TRUE, exact_summary_limit = 1e7) {
  
  ## Iterate over column names -------------------------------------------------
  col_list <- lapply(names(df), function(col) {
    vec <- df[[col]]
    
    ## basic descriptors -------------------------------------------------------
    is_num   <- is.numeric(vec)                   # numeric ? → mean/SD/hist
    cnt      <- length(vec)                       # total # of rows
    nulls    <- sum(is.na(vec))                    # NA count
    non_na_count <- cnt - nulls
    exact_summary <- cnt <= exact_summary_limit
    distinct_cnt <- if (exact_summary) data.table::uniqueN(vec, na.rm = TRUE) else NA_integer_
    
    ## numeric‑only statistics --------------------------------------------------
    has_numeric_values <- is_num && non_na_count > 0
    mn  <- if (has_numeric_values) min(vec, na.rm = TRUE) else NA_real_
    mx  <- if (has_numeric_values) max(vec, na.rm = TRUE) else NA_real_
    sov <- if (has_numeric_values) sum(vec, na.rm = na.rm) else NA_real_ # sum of values
    sosv <- if (has_numeric_values) sum((as.numeric(vec))^2, na.rm = na.rm) else NA_real_ # sum of squared values
    av  <- if (has_numeric_values) sov / non_na_count else NA_real_
    me <- if (has_numeric_values && exact_summary) median(vec, na.rm = na.rm) else NA_real_
    sdv <- if (has_numeric_values && non_na_count > 1) sqrt((sosv - (sov^2 / non_na_count)) / (non_na_count - 1)) else NA_real_
    n <- if (is_num) non_na_count else NA_real_# number of non NA values
    
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
