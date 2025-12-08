source("/src/parameters.R")
library(logger)
library(lidR)
source("/src/helper_functions.R")

######################################
## init
######################################

args <- parse_parameters()
output_folder <- args$output_dir
f <- args$dataset_path
m <- args$method

if(m == "single_file"){
  source("/src/run_single.R")
} else if(m == "collection"){
  source("/src/run_collection.R")
} else{
  stop("Method not recognized. Use 'single_file' or 'collection'")
}