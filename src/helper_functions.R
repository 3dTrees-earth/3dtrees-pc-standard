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
    geometry = sf::st_sfc(st_polygon(list(m)), crs = crs)
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
