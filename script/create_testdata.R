library(lidR)
library(rlas)
library(data.table)
library(sf)

# Create output folder
output_folder <- "faulty_las"
dir.create(output_folder, showWarnings = FALSE)

create_faulty_las <- function(filename, points, crs = NA, shift_bbox = FALSE, wrong_dtype = FALSE,
                              low_density = FALSE, small_extent = FALSE, bad_zrange = FALSE,
                              centroid_off_land = FALSE, high_point_format = FALSE) {
  # Basic grid of points
  X <- runif(points, 0, ifelse(small_extent, 0.5, 100))
  Y <- runif(points, 0, ifelse(small_extent, 0.5, 100))
  Z <- runif(points, ifelse(bad_zrange, -200, 0), ifelse(bad_zrange, 6000, 30))
  
  dt <- data.table(X = X, Y = Y, Z = Z)
  
  # Add an extra byte attribute with wrong R type
  if (wrong_dtype) {
    dt$Curvature <- as.character(sample(letters, points, replace = TRUE))
  } else {
    dt$Curvature <- as.numeric(runif(points))
  }
  
  # Lower density by creating large area
  if (low_density) {
    dt$X <- dt$X * 1000
    dt$Y <- dt$Y * 1000
  }
  
  # Build LAS object
  las <- LAS(dt)
  
  # Assign CRS
  if (!is.na(crs)) {
    st_crs(las) <- crs
  }
  
  # Shift bounding box to break bbox <-> point consistency
  if (shift_bbox) {
    header <- las@header
    header@PHB$`Min X` <- min(las@data$X) - 500
    header@PHB$`Max X` <- max(las@data$X) - 500
    header@PHB$`Min Y` <- min(las@data$Y) - 500
    header@PHB$`Max Y` <- max(las@data$Y) - 500
    las@header <- header
  }
  
  # Set Point Data Format
  if (high_point_format) {
    las@header$`Point Data Format ID` <- 8
    las@data$NIR <- sample(1000:5000, points, replace = TRUE)
  }
  
  # Save to file
  writeLAS(las, file.path(output_folder, filename))
}

# Faulty LAS examples
create_faulty_las("few_points.las", 50)  # < 100 points
create_faulty_las("no_crs.las", 1000, crs = NA)
create_faulty_las("bbox_mismatch.las", 1000, crs = 32632, shift_bbox = TRUE)
create_faulty_las("low_density.las", 1000, crs = 32632, low_density = TRUE)
create_faulty_las("small_extent.las", 1000, crs = 32632, small_extent = TRUE)
create_faulty_las("bad_zrange.las", 1000, crs = 32632, bad_zrange = TRUE)
create_faulty_las("wrong_extra_byte_type.las", 1000, crs = 32632, wrong_dtype = TRUE)
create_faulty_las("high_point_format.las", 1000, crs = 32632, high_point_format = TRUE)
create_faulty_las("non_faulty.las", 11000000, crs = 32632)
