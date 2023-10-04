#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/cohort_to_factors.R"))
source(here::here("src/shared/R-helpers/deseq2_analysis.R"))

suppressPackageStartupMessages({
  library("dplyr")
  library("readr")
})


settings <- get_settings(omic = "rnaseq", stage = "deg", dataset_name = NULL) #dataset_name="human"

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)


ensg_symbol <- readRDS(file.path(settings[["prefiltpca_dir"]], "ENSG_to_symbol.rds"))

# Helper functions for later using ensg_symbol
# character vector

# data frame with gene_id, will get gene_name column
df_geneid_to_gene_name <- function(x, add = FALSE) {
  if ("gene_name" %in% colnames(x)) {
    return(x)
  }
  out <- dplyr::left_join(x, ensg_symbol, by = "gene_id")
  if (add) {
    out <- dplyr::select(out, gene_id, gene_name, everything())
  } else {
    out <- dplyr::select(out, gene_name, everything())
  }
  out
}


cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = cols()
)


# Differential expression

cohort_for_deseq <- cohort_to_factors(cohort, settings[["cohort_factors"]])

rnaseq_intensities_mat_kept <- readRDS(file = file.path(settings[["prefiltpca_dir"]], "counts_raw_no_outliers_no_low_genes.rds"))

if (!all(colnames(rnaseq_intensities_mat_kept) == cohort$SampleID)) {
  stop("Mismatch of SampleIDs! This should not happen")
}

# Diff expr:

deseq2_analysis(
  count_data = rnaseq_intensities_mat_kept,
  col_data = cohort_for_deseq,
  deseq2_settings = settings[["deseq2"]],
  output_directory = settings[["output_directory"]],
  feature_column = "gene_id",
  results_callback = function(r) {
    r <- df_geneid_to_gene_name(r, add = TRUE)
    r
  }
)

