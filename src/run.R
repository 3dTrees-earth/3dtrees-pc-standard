source("parameters.R")
library(logger)

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
      stddev      = sdv
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

# setup logging i json format
logger::log_layout(logger::layout_json())

# args <- list(output_dir = getwd(), dataset_path = "C:/TLS/docker/input/Platane_pcrs.laz", min_density = 250)
args <- parse_parameters()
logger::log_threshold(INFO)
output_folder <- args$output_dir
f <- args$dataset_path
fname <- basename(f)
# strip off file ending
fname <- sub("\\.las$|\\.laz$", "", fname)
# send log to a json file named f.json in output dir
log_path <- file.path(output_folder, paste0(fname, ".json"))

logger::log_info("Logging to {log_path}")

if(file.exists(log_path)) {
  file.remove(log_path)
  logger::log_warn("Log {log_path} already existed and has been removed.")
}

logger::log_appender(appender_file(log_path))
logger::log_debug("Parameters parsed: {args}")

# load the necessary libraries
library(lidR)

# Try to read the file directly with lidR, bypassing extension check
logger::log_info("Attempting to read file directly: {f}")
logger::log_info("File size: {file.size(f)} bytes")

# load data 
# This is complete madness. It is necessary to bypass the extension check in lidR.
# This is a hack to get the tool to work in Galaxy.
logger::log_info("Loading data from {f}...")
las <- lidR::readLAS(f)

las_info_pre <- las_info(las)

if (is.null(las)) {
  stop("Failed to load LAS file with all methods")
}
logger::log_info("LAS file loaded successfully")
logger::log_info("Number of points: {nrow(las)}")
logger::log_info("Point format: {las@header$`Point Data Format ID`}")

# check the LAS object
check1 <- list(warnings = "Unable to run las_check!")
tryCatch(check1 <- lidR::las_check(las), "C++Error" = function(e){
  cat( sprintf( "C++ exception of class '%s' : %s\n", class(e)[1L], e$message  ) )
})


if(length(check1$errors) > 0) {
  logger::log_error("Errors found in the LAS object: {check1$errors}")
    #stop("Errors found in the LAS object")
}

if(length(check1$warnings) > 0) {
  logger::log_warn("Warnings found in the LAS object: {check1$warnings}")
}

if(length(check1$messages) > 0) {
  logger::log_warn("Messages found in the LAS object: {check1$messages}")
}


if(nrow(las) < 100){
  logger::log_error("Las file has < 100 points and will be dropped")
    stop("Las file has < 100 points and will be dropped")
} 

if(any(c("Invalid file: some points are outside the bounding box defined by the header","Invalid file: the bounding box defined by the header does not correspond to the points.") %in% check1$warnings)) {
    las <- lidR::las_update(las)
    logger::log_warn("bbox adapted since it did not match the points")
}

# check if all dimensions have an extend > 1m 
bbox <- sf::st_bbox(las)
extend_x <- bbox$xmax - bbox$xmin
extend_y <- bbox$ymax - bbox$ymin
extend_z <- max(las@data$Z, na.rm = TRUE) - min(las@data$Z, na.rm = TRUE)

logger::log_info("Extend in X: {extend_x}")
logger::log_info("Extend in Y: {extend_y}")
logger::log_info("Extend in Z: {extend_z}")

if(extend_x <= 1 || extend_y <= 1 || extend_z <= 1) {
  logger::log_error("Extend in one or more dimensions is less than or equal to 1m")
    stop("Extend in one or more dimensions is less than or equal to 1m")
}

# check if the Z coordinate range is plausible > -100 & < 5500 
if(any(extend_z < -100 || extend_z > 5500)) {
  logger::log_error("Z coordinate range is not plausible: < -100 or > 5500")
    stop("Z coordinate range is not plausible: < -100 or > 5500")
}

logger::log_info("Extend check compleate, start chull calculation.")

# calculate the point density
if(nrow(las) > 5e+08){ # point cloud > 500,000,000 points allocate gigantic amounts of memory
  idx <- sample(1:nrow(las), 5e+08)
  convex_hull <- las@data[idx,c("X","Y")][chull(las@data[idx,c("X","Y")]),] 
} else {
  convex_hull <- las@data[chull(las@data[,c("X","Y")]),c("X","Y")]
}

convex_hull <- list(as.matrix(rbind(convex_hull, convex_hull[1,]))) |> sf::st_polygon()
covered_area <- convex_hull |> sf::st_area()
point_density <- nrow(las@data) / covered_area
min_density_threshold <- args$min_density
logger::log_info("Point density: {point_density} points/m², minimum threshold: {min_density_threshold} points/m²")
if(point_density <= min_density_threshold) {
  logger::log_error("Point density is less than or equal to {min_density_threshold} points per square meter. Is this maybe ALS data or multiple separated plots in one file?")
    stop("Point density is less than or equal to ", min_density_threshold, " points per square meter. Is this maybe ALS data or multiple separated plots in one file?")
}

logger::log_info("Chull compleate, start attribute check.")

# check the extra byte arguments have the right data type from las@header@VLR$Extra_Bytes$`Extra Bytes Description`
# Type 1: char (1 byte) (not allowed in R)
# Type 2: unsigned char (1 byte) (not allowed in R)
# Type 3: short (2 bytes)
# Type 4: unsigned short (2 bytes)
# Type 5: long (4 bytes)
# Type 6: unsigned long (4 bytes)
# Type 7: long long (8 bytes)
# Type 8: unsigned long long (8 bytes)
# Type 9: float (4 byte IEEE floating point format)
# Type 10: double (8 byte IEEE floating point format)
# Type 11: string (a variable series of 1 byte characters, ASCII encoded, null-terminated (not allowed in R)
# mapping types to R types
types <- c("char", 
            "unsigned_char", 
            "integer", 
            "integer", 
            "integer", 
            "integer", 
            "integer", 
            "integer", 
            "numeric", 
            "numeric", 
            "string") 
#las@header@VLR$Extra_Bytes$`Extra Bytes Description`$Curvature$data_type <- 5
if(any(!(names(las@data) %in% lidR:::LASATTRIBUTES))){
    type_def <- sapply(las@header@VLR$Extra_Bytes$`Extra Bytes Description`, function(x) c(x["name"], x["data_type"])) |> t() |> data.frame()
    type_def$real_types <- sapply(las@data[,as.character(type_def$name), with = FALSE], class)
    
    if(!all(types[as.numeric(type_def$data_type)] == type_def$real_types)) {
      logger::log_warn("Extra byte arguments with the wrong data type detected")
    wrong_id <- which(types[as.numeric(type_def$data_type)] != type_def$real_types)
    for(i in wrong_id){
        #TODO: This is a hot-fix. Better check if a short or unsigned is also sufficient 
        las@header@VLR$Extra_Bytes$`Extra Bytes Description`[[i]]["data_type"] <- ifelse(type_def$real_types[i] == "numeric", 10,6)
    }
    }
}

# check for NA and static values in columns
removed <- 0
for(col in 1:ncol(las@data)){
  attribute <- colnames(las@data)[col-removed]
  
  # check for NA values 
  na_s <- which(is.na(las@data[,get(attribute)]))
  if(length(na_s) != 0) {
    if((length(na_s)/nrow(las@data)) > 0.01 & !(attribute %in% lidR:::LASATTRIBUTES)){
      las <- las |> remove_lasattribute(attribute)
      logger::log_warn("Las-attribute {attribute} contained NA values for more than 1% of points and had been removed.")
      removed <- removed + 1
      next
    } else if( !(attribute %in% lidR:::LASATTRIBUTES)) {
      las@data <- las@data[-na_s,]
      logger::log_warn("Las-attribute {attribute} contained NA values for less than 1% of points and those {length(na_s)} points have been removed.")
    } else {
      las@data[na_s, get(attribute)] <- 0
      logger::log_warn("Las-attribute {attribute} contained NA values for {length(na_s)} points and those have been set to 0.")
    }
  }
  
  # check for static values in columns
  if(length(unique(las@data[,get(attribute)])) == 1) {
    removed <- removed + 1
    if( !(attribute %in% lidR:::LASATTRIBUTES)) {
      las <- remove_lasattribute(las, attribute)
    } else {
      las@data <- las@data[,(attribute):=NULL]
    }
    logger::log_warn("Las-attribute {attribute} contained only one unique value and had been removed.")
  }
  
}


logger::log_info("Attributes compleate, start CRS check.")

# check if CRS is specified
has_crs <- !is.na(sf::st_crs(las))
log_info("CRS is not NA: {is.na(sf::st_crs(las))}")

# add to list 
if(!has_crs) {
  log_warn("CRS is NA")
  is_transformable <- FALSE
    # stop("CRS is NA")
} else {
    # is it possible to transform the CRS to 4978
    is_transformable <- tryCatch({
    random_sample <- rep(F, nrow(las@data))
    random_sample[sample(1:nrow(las@data), min(100, nrow(las@data)))] <- TRUE
    las_trans <- sf::st_transform(lidR::filter_poi(las, random_sample), 4978)
    TRUE
    }, error =  function(e) {
      log_error("Error in transforming CRS: {e$message}")
      FALSE
    })
    log_info("CRS is transformable to 4978: {is_transformable}")
    
    # convert convex_hull to EPSG:4326
    if(is_transformable){
      convex_hull_sf <- convex_hull |> sf::st_polygon() |> sf::st_sfc(crs = lidR::st_crs(las)) |> sf::st_transform(4326)
      wgs84_coordinates <- sf::st_coordinates(convex_hull_sf) |> as.data.frame()
      
      if(any(wgs84_coordinates$X >= 180 | wgs84_coordinates$X <= -180 | wgs84_coordinates$Y >= 90 | wgs84_coordinates$Y <= -90)){
        logger::log_warn("Implausible bbox detected.")
        is_transformable <- FALSE
      }
    }
}


# check if las Point Data Record Format is to high
if(las@header$`Point Data Format ID` > 7){
    if("NIR" %in% names(las@data)){
    las <- lidR::add_lasattribute(las, las@data$NIR, "NIRe", "the NIR channel changed to an extra byte attribute") #|> 
        #lidR::remove_lasattribute("NIR")
    las@data$NIR <- NULL
    }
    las@header$`Point Data Format ID` <- 7
    logger::log_warn("Las point data record format was >= 7 and has been downgraded. Full waveform might be lost. NIR might be renamed.")
}

# check if the bbox-centroid is on land area
if(is_transformable){
    #library(sf)
    bbox <- lidR::st_bbox(las)
    centroid <- sf::st_centroid(sf::st_as_sfc(bbox))
    
    # load land file
    land <- sf::st_read("geodata/ne_10m_land.shp", quiet = T)
    
    # check if the centroid is within the land area
    is_within <- length(sf::st_within(sf::st_transform(centroid,sf::st_crs(land)), land))
    logger::log_info("Centroid is within land area: {is_within > 0}")
    
    # second try with small islands
    if(is_within == 0) {
    land <- sf::st_read("geodata/ne_10m_minor_islands.shp", quiet = T)
    is_within <- length(sf::st_within(sf::st_transform(centroid,sf::st_crs(land)), land))
    logger::log_info("Centroid is within minor island area: {is_within > 0}")
    }
    
    if(is_within == 0) {
      logger::log_error("Centroid is not within land area")
      is_transformable <- FALSE
      #stop("Centroid is not within land area")
    }
} 
if(!is_transformable){
    logger::log_error("CRS is not transformable to epsg:4978")
    logger::log_warn("Reoffset dataset to 0,0,0, fake-CRS will be set to EPSG:4978")
    las@data$X <- las@data$X - min(las@data$X)
    las@data$Y <- las@data$Y - min(las@data$Y)
    las@data$Z <- las@data$Z - min(las@data$Z)
    las <- las |> lidR::las_reoffset()
    lidR::st_crs(las) <- 4978
    # stop("CRS is not transformable to epsg:3857")
}

logger::log_info("CRS compleate, write outputs.")


# summarystats <- function(x) c(min = min(x), max = max(x), mean = mean(x), median = median(x), sd = sd(x), q = quantile(x, c(0.05, 0.95)), n_unique = length(unique(x)),unique_point_ratio = length(unique(x)) / length(x), type = typeof(x))
# col_info <- list()

has_RGB <- all(c("R","G","B") %in% colnames(las@data))
# col_info$summary_labels <- c("min","max","mean","median","sd","qantile5percent","qantile95percent", "n_unique", "unique_point_ratio", "type")
# col_info$attribute_summary <- lapply(las@data, summarystats)
# 
# 
# con <- file(log_path, "a+")
#   jsonlite::write_json(col_info, con)
# close(con)
logger::log_info("Point cloud contains RGB values: {has_RGB}")

# logger::log_info("Color information compleate, write outputs.")

las_info_post <- las_info(las)

# change directory to output folder
setwd(output_folder)

# write the las file to the output folder
output_file <- "pc_standardized.laz"
writeLAS(las, output_file)
if(is_transformable){
  
  sf::st_write(convex_hull_sf, "convex_hull_wgs84.GeoJSON", append=FALSE)
  logger::log_info("CRS is valid and transformable and the convex hull was written to convex_hull_wgs84.GeoJSON")
}
logger::log_info("File {fname} is valid and has been written to {output_file}.")
logger::log_info("Writing Lasinfo for pointcloud before and after standardization.")
las_info_pre$standardized <- FALSE
las_info_post$standardized <- TRUE
con <- file(log_path, "a+")
  jsonlite::write_json(las_info_pre, con)
  jsonlite::write_json(las_info_post, con)
close(con)
