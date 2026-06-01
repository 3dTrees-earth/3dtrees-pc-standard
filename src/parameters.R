parse_approx_location <- function(value) {
    if (is.null(value)) {
        return(NULL)
    }

    trimmed <- trimws(as.character(value)[1])
    if (!nzchar(trimmed)) {
        return(NULL)
    }

    number_pattern <- "[-+]?(?:[0-9]*\\.[0-9]+|[0-9]+\\.?[0-9]*)(?:[eE][-+]?[0-9]+)?"
    matches <- regmatches(trimmed, gregexpr(number_pattern, trimmed, perl = TRUE))[[1]]

    if (length(matches) != 2) {
        stop("Approximate location must contain exactly two numbers in lon, lat order, for example: (7.85, 47.99)")
    }

    lon <- suppressWarnings(as.numeric(matches[[1]]))
    lat <- suppressWarnings(as.numeric(matches[[2]]))

    if (!is.finite(lon) || lon < -180 || lon > 180) {
        stop("Approximate location longitude must be between -180 and 180.")
    }
    if (!is.finite(lat) || lat < -90 || lat > 90) {
        stop("Approximate location latitude must be between -90 and 90.")
    }

    list(lon = lon, lat = lat)
}

validate_approx_coordinate_pair <- function(lon, lat) {
    lon_missing <- is.null(lon) || length(lon) == 0 || is.na(lon)
    lat_missing <- is.null(lat) || length(lat) == 0 || is.na(lat)

    if (xor(lon_missing, lat_missing)) {
        stop("Approximate longitude and latitude must be provided together.")
    }

    if (!lon_missing && !lat_missing) {
        if (!is.finite(lon) || lon < -180 || lon > 180) {
            stop("Approximate longitude must be between -180 and 180.")
        }
        if (!is.finite(lat) || lat < -90 || lat > 90) {
            stop("Approximate latitude must be between -90 and 90.")
        }
    }
}

parse_parameters <- function() {
    library(argparse)

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
        "--approx-location", "--approx_location",
        type = "character",
        required = FALSE,
        help = "Approximate collection WGS84 point as a lon, lat tuple, for example: (7.85, 47.99). Overrides the calculated collection centroid."
    )

    parser$add_argument(
        "--approx-lon", "--approx_lon",
        type = "double",
        required = FALSE,
        help = "Deprecated. Approximate collection longitude in EPSG:4326. Use --approx-location instead."
    )

    parser$add_argument(
        "--approx-lat", "--approx_lat",
        type = "double",
        required = FALSE,
        help = "Deprecated. Approximate collection latitude in EPSG:4326. Use --approx-location instead."
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

    approx_location <- parse_approx_location(args$approx_location)
    if (!is.null(approx_location)) {
        if (!is.null(args$approx_lon) || !is.null(args$approx_lat)) {
            stop("Use either --approx-location or --approx-lon/--approx-lat, not both.")
        }
        args$approx_lon <- approx_location$lon
        args$approx_lat <- approx_location$lat
    }

    validate_approx_coordinate_pair(args$approx_lon, args$approx_lat)

    return(args)
}
