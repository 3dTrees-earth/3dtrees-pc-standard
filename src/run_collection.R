######################################
## collection check
######################################

files <- list.files(
  f,
  pattern = "\\.la[sz]$",
  full.names = TRUE,
  ignore.case = TRUE
)

files_info <- list()

i <- 1
for(file in files){
  file_info <- list()
  file_info$id <- i
  file_info$f_name <- basename(file)
  las <- readLAS(file)
  # get the crs
  file_info$crs <- st_crs(las)
  # get column names of first tile 
  file_info$attributes <- colnames(las@data)
  # get collum types
  file_info$attr_types <- sapply(las@data, class)
  
  # calculate convex hull
  if(nrow(las) > 5e+08){ # point cloud > 500,000,000 points allocate gigantic amounts of memory
    idx <- sample(1:nrow(las), 5e+08)
    file_info$concave_hull <- las@data[idx,c("X","Y")][chull(las@data[idx,c("X","Y")]),] 
  } else {
    file_info$concave_hull <- las@data[chull(las@data[,c("X","Y")]),c("X","Y")]
  }
  
  # number of points
  file_info$n_points <- nrow(las@data)
  
  # stats on the attributes including sum and squared sum for SD and mean calculation over all tiles
  file_info$attribute_stats <- attribute_stats(las@data)
  
  # write to list with all tile infos
  files_info[[file_info$f_name]] <- file_info
  
  # increase counter
  i <- i+1
}

######################################
## collection‑level checks
######################################

## 1) CRS homogeneity ---------------------------------------------------------
crs_list <- lapply(files_info, `[[`, "crs")
crs_epsg <- sapply(crs_list, function(x) x$epsg)
unique_epsg <- unique(crs_epsg)
homogeneous_crs <- length(unique_epsg) == 1

## 2) Attribute names homogeneity --------------------------------------------
attr_names_list <- lapply(files_info, `[[`, "attributes")

# Reference attribute names (from first tile)
ref_attr_names  <- attr_names_list[[1]]

# Check whether each tile has exactly the same attribute vector as the first
same_attr_names <- vapply(
  attr_names_list,
  function(x) identical(x, ref_attr_names),
  logical(1)
)
homogeneous_attr_names <- all(same_attr_names)

# All attribute names across all tiles (union)
all_attr_names <- sort(unique(unlist(attr_names_list)))

# Attributes present in *all* tiles (intersection)
common_attr_names <- Reduce(intersect, attr_names_list)

# Attributes that are NOT present in all tiles ("non overlapping")
non_overlapping_attr_names <- setdiff(all_attr_names, common_attr_names)


# attribute_stats stored per file in files_info[[...]]$attribute_stats
attr_stat_list <- lapply(files_info, `[[`, "attribute_stats")

only_na_attr <- c()

for (attr_nm in all_attr_names) {
  # collect non‑NA counts n for this attribute across tiles
  n_vec <- numeric(length(attr_stat_list))
  
  for (k in seq_along(attr_stat_list)) {
    tile_stats <- attr_stat_list[[k]]
    idx <- which(vapply(tile_stats, `[[`, character(1), "name") == attr_nm)
    if (length(idx) == 1) {
      n_vec[k] <- tile_stats[[idx]]$n   # number of non‑NA values in this tile
    } else {
      n_vec[k] <- 0
    }
  }
  
  # if total non‑NA count is zero, attribute is all‑NA across all tiles
  if (all(is.na(n_vec))) {
    only_na_attr <- c(only_na_attr, attr_nm)
  }
}

zero_mean_sd_attr <- c()

for (attr_nm in all_attr_names) {
  all_zero <- TRUE
  
  for (k in seq_along(attr_stat_list)) {
    tile_stats <- attr_stat_list[[k]]
    idx <- which(vapply(tile_stats, `[[`, character(1), "name") == attr_nm)
    
    if (length(idx) == 1) {
      st <- tile_stats[[idx]]
      
      # only consider tiles where this attribute has at least one non‑NA value
      if (isTRUE(st$n > 0)) {
        # here we check for exactly 0; if you want a tolerance use abs(...) < eps
        if (!(st$average == 0 && st$stddev == 0)) {
          all_zero <- FALSE
          break
        }
      }
    }
  }
  
  if (all_zero) {
    zero_mean_sd_attr <- c(zero_mean_sd_attr, attr_nm)
  }
}

non_na_attr <- setdiff(all_attr_names, only_na_attr)
non_na_0_attr <- setdiff(non_na_attr, zero_mean_sd_attr)

## 3) Attribute types homogeneity --------------------------------------------
attr_types_list <- lapply(files_info, `[[`, "attr_types")
ref_attr_types  <- attr_types_list[[1]]
same_attr_types <- vapply(
  attr_types_list,
  function(x) identical(x[names(ref_attr_types)], ref_attr_types),
  logical(1)
)
homogeneous_attr_types <- all(same_attr_types)

## 4) Polygons from convex hulls + overlaps / disjoint ------------------------

# per‑tile POLYGONs (using helper)
tile_polys <- lapply(
  files_info,
  function(info) {
    chull_to_polygon(
      coords = info$concave_hull,
      crs    = info$crs
    )
  }
)

# bind to one sf with a "tile" attribute
tile_polys <- Map(
  function(poly, nm) {
    poly$tile <- nm
    poly
  },
  tile_polys,
  names(files_info)
)
tiles_sf <- do.call(rbind, tile_polys)

#TODO: add 1 m buffers as safety margins

# overlaps (true overlap, not just touching)
ov_mat <- sf::st_overlaps(sf::st_buffer(tiles_sf,-1), tiles_sf, sparse = FALSE)
diag(ov_mat) <- FALSE
any_overlap <- any(ov_mat)

# disjoint?
disj_mat <- sf::st_disjoint(sf::st_buffer(tiles_sf,100), tiles_sf, sparse = FALSE)
diag(disj_mat) <- TRUE
all_disjoint <- all(disj_mat)
any_disjoint <- any(!apply(disj_mat,1,any))

## 5) Multipolygon WKT (concave hull map) -------------------------------------

mp_sf  <- polygons_to_multipolygon(tile_polys, files = names(files_info)) |> sf::st_geometry()
sf::st_crs(mp_sf) <- crs_list[[1]]
if(is_crs_defined(crs_list[[1]]) & is_crs_transformable(crs_list[[1]])){
  mp_sf <- sf::st_transform(mp_sf, 4326)
}

mp_wkt <- sf::st_as_text(mp_sf[[1]])

## 6) Global attribute statistics over all tiles -------------------------------

attr_stat_list <- lapply(files_info, `[[`, "attribute_stats")

global_attr_stats <- lapply(non_na_0_attr, function(attr_nm) {
  n_vec    <- numeric(length(attr_stat_list))
  sov_vec  <- numeric(length(attr_stat_list))
  sosv_vec <- numeric(length(attr_stat_list))
  
  for (k in seq_along(attr_stat_list)) {
    tile_stats <- attr_stat_list[[k]]
    idx <- which(vapply(tile_stats, `[[`, character(1), "name") == attr_nm)
    if (length(idx) == 1) {
      n_vec[k]    <- tile_stats[[idx]]$n
      sov_vec[k]  <- tile_stats[[idx]]$sumofvals
      sosv_vec[k] <- tile_stats[[idx]]$sumofsqvals
    } else {
      n_vec[k]    <- 0
      sov_vec[k]  <- 0
      sosv_vec[k] <- 0
    }
  }
  
  mtv <- multi_tile_variance(n_vec, sov_vec, sosv_vec)
  
  list(
    name     = attr_nm,
    n_total  = sum(n_vec),
    mean     = mtv$mean,
    variance = mtv$variance,
    stddev   = mtv$stddev
  )
})

######################################
## prepare JSON structure
######################################

# make CRS JSON‑friendly (convert st_crs objects to simple lists)
files_info_json <- lapply(files_info, function(info) {
  info$crs <- as.list(info$crs)
  attr(info$crs, "class") <- NULL
  info
})

collection_checks <- list(
  n_tiles                      = length(files_info),
  homogeneous_crs              = homogeneous_crs,
  first_crs                    = files_info_json[[1]]$crs,
  crs_epsg_values              = unique_epsg,
  homogeneous_attribute_names  = homogeneous_attr_names,
  reference_attribute_names    = ref_attr_names,
  homogeneous_attribute_types  = homogeneous_attr_types,
  reference_attribute_types    = ref_attr_types,
  tiles_any_overlap            = any_overlap,
  tiles_all_disjoint           = all_disjoint,
  global_attribute_stats       = global_attr_stats,
  multipolygon_wkt             = mp_wkt,
  all_attribute_names          = all_attr_names,
  common_attribute_names       = common_attr_names,
  non_overlapping_attributes   = non_overlapping_attr_names,
  removeable_attributes        = c(only_na_attr, zero_mean_sd_attr)
)

result <- list(
  collection = collection_checks,
  files      = files_info_json
)

out_file <- file.path(output_folder, "collection_summary.json")
jsonlite::write_json(result, out_file, pretty = TRUE, auto_unbox = TRUE, na = "null")
