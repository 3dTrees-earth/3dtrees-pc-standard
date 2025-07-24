# validate_las.R

input_folder <- "/input/"
output_folder <- "/output/"

# get all las and laz files
files <- list.files(input_folder, pattern = "\\.las$|\\.laz$", full.names = TRUE)

# load the necessary libraries 
library(lidR)

# loop over files
for(f in files){
  
  fname <- basename(f)
  # strip off file ending
  fname <- sub("\\.las$|\\.laz$", "", fname)
  
  # load data
  las <- lidR::readLAS(f)
  
  # create empty list for warnings and errors
  warns <- list()
  errors <- list()
  
  # check the LAS object
  check1 <- lidR::las_check(las)
  
  if(length(check1$errors) > 0) {
    errors <- c(errors, check1$errors)
  }
  
  if(nrow(las) < 100){
    errors <- c(errors, "Las file has < 100 points and will be dropped")
  } else {
    
    if(any(c("Invalid file: some points are outside the bounding box defined by the header","Invalid file: the bounding box defined by the header does not correspond to the points.") %in% check1$warnings)) {
      las <- lidR::las_update(las)
      warns <- c(warns, "bbox adapted since it did not match the points")
    }
    
    # check if all dimensions have an extend > 1m 
    
    bbox <- sf::st_bbox(las)
    extend_x <- bbox$xmax - bbox$xmin
    extend_y <- bbox$ymax - bbox$ymin
    extend_z <- max(las@data$Z, na.rm = TRUE) - min(las@data$Z, na.rm = TRUE)
    
    print(paste("Extend in X:", extend_x))
    print(paste("Extend in Y:", extend_y))
    print(paste("Extend in Z:", extend_z))
    
    if(extend_x <= 1 || extend_y <= 1 || extend_z <= 1) {
      errors <- c(errors, "Extend in one or more dimensions is less than or equal to 1m")
    }
    
    # check if the Z coordinate range is plausible > -100 & < 5500 
    if(any(extend_z < -100 || extend_z > 5500)) {
      errors <- c(errors, "Z coordinate range is not plausible: < -100 or > 5500")
    }
    
    # calculate the point density
    convex_hull <- las@data[chull(las@data[,c("X","Y")]),c("X","Y")]
    covered_area <- list(as.matrix(rbind(convex_hull, convex_hull[1,]))) |> sf::st_polygon() |> sf::st_area()
    point_density <- nrow(las@data) / covered_area
    if(point_density <= 1000) {
      errors <- c(errors, "Point density is less than or equal to 1000 points per square meter. Is this maybe ALS data or multiple separated plots in one file?")
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
        warns <- c(warns, "Extra byte arguments with the wrong data type detected")
        wrong_id <- which(types[as.numeric(type_def$data_type)] != type_def$real_types)
        for(i in wrong_id){
          #TODO: This is a hot-fix. Better check if a short or unsigned is also sufficient 
          las@header@VLR$Extra_Bytes$`Extra Bytes Description`[[i]]["data_type"] <- ifelse(type_def$real_types[i] == "numeric", 10,6)
        }
      }
    }
    
    # check if CRS is specified
    has_crs <- !is.na(sf::st_crs(las))
    print(paste("CRS is not NA:", is.na(sf::st_crs(las))))
    
    # add to list 
    if(!has_crs) {
      errors <- c(errors, "CRS is NA")
    } else {
      # is it possible to transform the CRS to 3857
      is_transformable <- tryCatch({
        random_sample <- rep(F, nrow(las@data))
        random_sample[sample(1:nrow(las@data), min(100, nrow(las@data)))] <- TRUE
        las_trans <- sf::st_transform(lidR::filter_poi(las, random_sample), 3857)
        TRUE
      }, error =  function(e) {
        message("Error in transforming CRS: ", e$message)
        FALSE
      })
      print(paste("CRS is transformable to 3857:", is_transformable))
      
      # add to list
      if(!is_transformable) {
        errors <- c(errors, "CRS is not transformable to 3857")
      }
    }
    
    
    # check if las Point Data Record Format is to high
    if(las@header$`Point Data Format ID` > 7){
      if("NIR" %in% names(las@data)){
        las <- lidR::add_lasattribute(las, las@data$NIR, "NIRe", "the NIR channel changed to an extra byte attribute") |> 
          lidR::remove_lasattribute("NIR")
      }
      las@header$`Point Data Format ID` <- 7
      warns <- c(warns, "Las point data record format was >= 7 and has been downgraded. Full waveform might be lost. NIR might be renamed.")
    }
    
    # check if the bbox-centroid is on land area
    if(has_crs && is_transformable){
      library(sf)
      bbox <- lidR::st_bbox(las)
      centroid <- sf::st_centroid(sf::st_as_sfc(bbox))
      
      # load land file
      land <- sf::st_read("geodata/ne_10m_land.shp", quiet = T)
      
      # check if the centroid is within the land area
      is_within <- length(sf::st_within(sf::st_transform(centroid,sf::st_crs(land)), land))
      print(paste("Centroid is within land area:", is_within > 0))
      
      # second try with small islands
      if(is_within == 0) {
        land <- sf::st_read("geodata/ne_10m_minor_islands.shp", quiet = T)
        is_within <- length(sf::st_within(sf::st_transform(centroid,sf::st_crs(land)), land))
        print(paste("Centroid is within minor island area:", is_within > 0))
      }
      
      if(is_within == 0) {
        errors <- c(errors, "Centroid is not within land area")
      }
    }
  }
  # write warns and errors to one json file
  json <- jsonlite::toJSON(list(list(warnings = warns), errors = errors))
  jsonlite::write_json(json, paste0(output_folder, fname, "_protocol.json") )
  
  # write the las file to the output folder
  if(length(errors) == 0) {
    writeLAS(las, paste0(output_folder, fname,".laz"))
    message(paste("File", fname, "is valid and has been written to the output folder."))
  } else {
    message(paste("File", fname, "is not valid and has not been written to the output folder."))
  }
}
  
