#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/cohort_to_factors.R"))
source(here::here("src/shared/R-helpers/deseq2_analysis.R"))

suppressPackageStartupMessages({
  library("readr")
})

settings <- get_settings(omic = "srna", stage = "deg", dataset_name = NULL) # dataset_name="c9")

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)


cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = cols()
)


cohort_for_deseq <- cohort_to_factors(cohort, settings[["cohort_factors"]])


mature_to_keep <- readRDS(file.path(settings[["prefilt_dir"]], "to_keep_mirna_mature.rds"))
mature_counts_mat <- readRDS(file = file.path(settings[["prefilt_dir"]], "mature_counts_mat.rds"))
if (!all(colnames(mature_counts_mat) == cohort$SampleID)) {
  stop("Mismatch of SampleIDs! This should not happen")
}

hairpin_to_keep <- readRDS(file.path(settings[["prefilt_dir"]], "to_keep_mirna_hairpin.rds"))
hairpin_counts_mat <- readRDS(file = file.path(settings[["prefilt_dir"]], "hairpin_counts_mat.rds"))
if (!all(colnames(hairpin_counts_mat) == cohort$SampleID)) {
  stop("Mismatch of SampleIDs! This should not happen")
}


# Diff expr:

deseq2_analysis(
  count_data = mature_counts_mat[mature_to_keep,],
  col_data = cohort_for_deseq,
  deseq2_settings = settings[["deseq2"]],
  output_directory = file.path(settings[["output_directory"]], "mature"),
  feature_column = "mirnaID",
  results_callback = NULL
)


deseq2_analysis(
  count_data = hairpin_counts_mat[hairpin_to_keep,],
  col_data = cohort_for_deseq,
  deseq2_settings = settings[["deseq2"]],
  output_directory = file.path(settings[["output_directory"]], "hairpin"),
  feature_column = "mirnaID",
  results_callback = NULL
)