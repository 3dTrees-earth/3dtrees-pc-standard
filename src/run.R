source("parameters.R")
library(logger)

# setup logging i json format
logger::log_layout(logger::layout_json())

args <- parse_parameters()
log_threshold(DEBUG)
log_debug("Parameters parsed: {args}")

# load the necessary libraries
library(lidR)

# Try to read the file directly with lidR, bypassing extension check
log_info("Attempting to read file directly: {args$dataset_path}")
log_info("File size: {file.size(args$dataset_path)} bytes")

# copy from validate_las.R
# validate_las.R

output_folder <- args$output_dir
f <- args$dataset_path

# send log to stderr instead of a file to avoid permission issues
logger::log_appender(appender_stderr)

  
fname <- basename(f)
# strip off file ending
fname <- sub("\\.las$|\\.laz$", "", fname)

# load data 
# This is complete madness. It is necessary to bypass the extension check in lidR.
# This is a hack to get the tool to work in Galaxy.
log_info("Loading data from {args$dataset_path}...")
las <- lidR::readLAS(args$dataset_path)

if (is.null(las)) {
  stop("Failed to load LAS file with all methods")
}
log_info("LAS file loaded successfully")
log_info("Number of points: {nrow(las)}")
log_info("Point format: {las@header$`Point Data Format ID`}")

# check the LAS object
check1 <- lidR::las_check(las)

if(length(check1$errors) > 0) {
    log_error("Errors found in the LAS object: {check1$errors}")
    #stop("Errors found in the LAS object")
}

if(length(check1$warnings) > 0) {
    log_warn("Warnings found in the LAS object: {check1$warnings}")
}

if(nrow(las) < 100){
    log_error("Las file has < 100 points and will be dropped")
    stop("Las file has < 100 points and will be dropped")
} 

if(any(c("Invalid file: some points are outside the bounding box defined by the header","Invalid file: the bounding box defined by the header does not correspond to the points.") %in% check1$warnings)) {
    las <- lidR::las_update(las)
    log_warn("bbox adapted since it did not match the points")
}

# check if all dimensions have an extend > 1m 
bbox <- sf::st_bbox(las)
extend_x <- bbox$xmax - bbox$xmin
extend_y <- bbox$ymax - bbox$ymin
extend_z <- max(las@data$Z, na.rm = TRUE) - min(las@data$Z, na.rm = TRUE)

log_info("Extend in X: {extend_x}")
log_info("Extend in Y: {extend_y}")
log_info("Extend in Z: {extend_z}")

if(extend_x <= 1 || extend_y <= 1 || extend_z <= 1) {
    log_error("Extend in one or more dimensions is less than or equal to 1m")
    stop("Extend in one or more dimensions is less than or equal to 1m")
}

# check if the Z coordinate range is plausible > -100 & < 5500 
if(any(extend_z < -100 || extend_z > 5500)) {
    log_error("Z coordinate range is not plausible: < -100 or > 5500")
    stop("Z coordinate range is not plausible: < -100 or > 5500")
}

# calculate the point density
convex_hull <- las@data[chull(las@data[,c("X","Y")]),c("X","Y")]
covered_area <- list(as.matrix(rbind(convex_hull, convex_hull[1,]))) |> sf::st_polygon() |> sf::st_area()
point_density <- nrow(las@data) / covered_area
if(point_density <= 500) {
    log_error("Point density is less than or equal to 500 points per square meter. Is this maybe ALS data or multiple separated plots in one file?")
    stop("Point density is less than or equal to 500 points per square meter. Is this maybe ALS data or multiple separated plots in one file?")
}

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
    log_warn("Extra byte arguments with the wrong data type detected")
    wrong_id <- which(types[as.numeric(type_def$data_type)] != type_def$real_types)
    for(i in wrong_id){
        #TODO: This is a hot-fix. Better check if a short or unsigned is also sufficient 
        las@header@VLR$Extra_Bytes$`Extra Bytes Description`[[i]]["data_type"] <- ifelse(type_def$real_types[i] == "numeric", 10,6)
    }
    }
}

# check if CRS is specified
has_crs <- !is.na(sf::st_crs(las))
log_info("CRS is not NA: {is.na(sf::st_crs(las))}")

# add to list 
if(!has_crs) {
    log_warn("CRS is NA")
    # stop("CRS is NA")
} else {
    # is it possible to transform the CRS to 3857
    is_transformable <- tryCatch({
    random_sample <- rep(F, nrow(las@data))
    random_sample[sample(1:nrow(las@data), min(100, nrow(las@data)))] <- TRUE
    las_trans <- sf::st_transform(lidR::filter_poi(las, random_sample), 3857)
    TRUE
    }, error =  function(e) {
    log_error("Error in transforming CRS: {e$message}")
    FALSE
    })
    log_info("CRS is transformable to 3857: {is_transformable}")
    
    # add to list
    if(!is_transformable) {
    log_error("CRS is not transformable to epsg:3857")
    # stop("CRS is not transformable to epsg:3857")
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
    log_warn("Las point data record format was >= 7 and has been downgraded. Full waveform might be lost. NIR might be renamed.")
}

# check if the bbox-centroid is on land area
if(has_crs && is_transformable){
    #library(sf)
    bbox <- lidR::st_bbox(las)
    centroid <- sf::st_centroid(sf::st_as_sfc(bbox))
    
    # load land file
    land <- sf::st_read("geodata/ne_10m_land.shp", quiet = T)
    
    # check if the centroid is within the land area
    is_within <- length(sf::st_within(sf::st_transform(centroid,sf::st_crs(land)), land))
    log_info("Centroid is within land area: {is_within > 0}")
    
    # second try with small islands
    if(is_within == 0) {
    land <- sf::st_read("geodata/ne_10m_minor_islands.shp", quiet = T)
    is_within <- length(sf::st_within(sf::st_transform(centroid,sf::st_crs(land)), land))
    log_info("Centroid is within minor island area: {is_within > 0}")
    }
    
    if(is_within == 0) {
        log_error("Centroid is not within land area")
        stop("Centroid is not within land area")
    }
}

# change directory to output folder
setwd(output_folder)

# write the las file to the output folder
output_file <- "pc_standardized.laz"
writeLAS(las, output_file)
log_info("File {fname} is valid and has been written to {output_file}.")

  
