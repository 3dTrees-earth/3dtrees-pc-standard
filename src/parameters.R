library(argparse)

parse_parameters <- function() {
    # Setup the parser
    parser <- ArgumentParser(description = "3DTrees LAS/LAZ validation and standardization tool")

    parser$add_argument(
        "--dataset-path", "--dataset_path",
        type = "character",
        required = TRUE,
        default = "/in/",
        help = "Input dataset path of (LAS/LAZ) files folder"
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