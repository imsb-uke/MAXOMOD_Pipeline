suppressPackageStartupMessages({
  requireNamespace("argparse")
})

get_args <- function() {

  # create parser object
  parser <- argparse::ArgumentParser()

  parser$add_argument(
    "--dataset-name",
    type = "character",
    default = NULL,
    help = "dataset name we apply the script to (e.g. sod1)"
  )

  args <- parser$parse_args()
  args
}
