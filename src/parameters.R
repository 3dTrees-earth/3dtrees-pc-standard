library(argparse)

parse_parameters <- function() {
    # Setup the parser
    parser <- ArgumentParser(description = "3DTrees LAS/LAZ validation and standardization tool")
    
    parser$add_argument(
        "--dataset-path", "--dataset_path",
        type = "character",
        required = TRUE,
        default = "/in/",
        help = "Input dataset path of (LAS/LAZ) files or folders"
    )

    parser$add_argument(
        "--output-dir", "--output_dir",
        type = "character",
        required = TRUE,
        default = "/out",
        help = "Output directory"
    )

    parser$add_argument(
        "--min-density", "--min_density",
        type = "integer",
        required = FALSE,
        default = 10,
        help = "Minimum point density threshold (points per square meter)"
    )
    
    parser$add_argument(
        "--removeable_attributes", "--removeable_attributes",
        type = "character",
        required = FALSE,
        default = '',
        nargs='+',
        help = "Argument to remove collumns, usually NA or empty ones, for single file checks."
    )
    
    parser$add_argument(
        "--method",
        type = "character",
        default = "single_file",
        help = "Is the input a single las/laz file or a colelction. Collection will check the collection while single files get standardized."
    )

    parser$add_argument(
        "--approx-lon", "--approx_lon",
        type = "double",
        required = FALSE,
        help = "Approximate collection longitude in EPSG:4326. When paired with --approx-lat, overrides the calculated collection centroid."
    )

    parser$add_argument(
        "--approx-lat", "--approx_lat",
        type = "double",
        required = FALSE,
        help = "Approximate collection latitude in EPSG:4326. When paired with --approx-lon, overrides the calculated collection centroid."
    )

    # Parse the arguments
    args <- parser$parse_args()

    # Validate the arguments
    if (!file.exists(args$dataset_path)) {
        stop("Input dataset path does not exist: ", args$dataset_path)
    }

    if (!dir.exists(args$output_dir)) {
        stop("Output directory does not exist: ", args$output_dir)
    }

    return(args)
}
