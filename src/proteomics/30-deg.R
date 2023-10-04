#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/cohort_to_factors.R"))
source(here::here("src/shared/R-helpers/limma_analysis.R"))

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("readr")
  library("tibble")
})

settings <- get_settings(omic = "proteomics", stage = "deg", dataset_name = NULL) # dataset_name = "sod1"

cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = readr::cols()
)

intensity_mat_for_da <- readRDS(settings[["intensity_mat_for_limma"]])

if (!all(rownames(intensity_mat_for_da) == cohort$SampleID)) {
  stop("Error")
}

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

cohort_for_limma <- cohort_to_factors(cohort, settings[["cohort_factors"]])

de_results <- limma_analysis(
  x = intensity_mat_for_da,
  col_data = cohort_for_limma,
  limma_settings = settings[["limma"]],
  output_directory = settings[["output_directory"]],
  results_callback = NULL
)

saveRDS(de_results, file = file.path(settings[["output_directory"]], "de_results.rds"))
