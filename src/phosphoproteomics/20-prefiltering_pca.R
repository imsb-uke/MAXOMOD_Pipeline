#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/pca_fit.R"))
source(here::here("src/shared/R-helpers/pc_scores_loadings_plot.R"))

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("readr")
  library("tibble")
  library("ggplot2")
  library("missForest")
  requireNamespace("forcats")
  library("preprocessCore")
})


settings <- get_settings(omic = "phosphoproteomics", stage = "prefiltering_pca", dataset_name = NULL) # dataset_name = "sod1"

cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = NULL
)


disease_grp <- setdiff(cohort$Condition, "ctrl")
colours_sex <- c("female" = "#df65b0", "male" = "#2c7fb8")
colours_condition <- c("ctrl" = "#fe9929", "als" = "#31a354")
names(colours_condition) <- c("ctrl", disease_grp)


dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)


intensity <- readr::read_csv(file.path(settings[["input_intensity"]]), col_types = NULL)
intensity <- column_to_rownames(intensity, var = "SampleID")
intensity_mat <- as.matrix(intensity)
rownames(intensity_mat) <- rownames(intensity)

if (!all(rowSums(is.na(intensity_mat)) == 0)) {
  stop("Missing values detected. Unexpected.")
}

intensity_mat <- intensity_mat[cohort$SampleID, ]

saveRDS(intensity_mat, file = file.path(settings[["output_directory"]], "intensity_mat_filtered_imputed_log2transf_norm.rds"))


# PCA --------------

pcadir <- file.path(settings[["output_directory"]], "pca")
dir.create(pcadir, recursive = TRUE, showWarnings = FALSE)

message("Creating PCA")
pca_model_info <- pca_fit(intensity_mat, ntop = 500, cent = TRUE, scal = FALSE, scores_metadata = cohort)
saveRDS(pca_model_info, file = file.path(pcadir, "pca_model_info.rds"))

message("PCA plots")

for (pca_colour_by in settings[["pca_colour_by"]]) {
  if (pca_colour_by == "Condition") {
    color_values <- colours_condition
  } else {
    color_values <- NULL
  }
  pc_12_plt <- pc_scores_loadings_plot(
    pca_model_info,
    pcx = 1, pcy = 2,
    color = pca_colour_by,
    color_values = color_values,
    loading_label = "Protein"
  )
  saveRDS(pc_12_plt, file = file.path(pcadir, sprintf("pca_12_%s.rds", pca_colour_by)))
  ggsave(
    filename = file.path(pcadir, sprintf("pca_12_%s.pdf", pca_colour_by)),
    plot = pc_12_plt,
    dpi = 300, width = unit(10, "in")
  ) # , height = unit(5, "in"))

  pc_34_plt <- pc_scores_loadings_plot(
    pca_model_info,
    pcx = 3, pcy = 4,
    color = pca_colour_by, color_values = color_values,
    loading_label = "Protein"
  )

  saveRDS(pc_34_plt, file = file.path(pcadir, sprintf("pca_34_%s.rds", pca_colour_by)))
  ggsave(
    filename = file.path(pcadir, sprintf("pca_34_%s.pdf", pca_colour_by)),
    plot = pc_34_plt,
    dpi = 300, width = unit(10, "in") # , height = unit(5, "in")
  )
}

# PCA - all vars --------------
message("Creating PCA")
pca_model_info <- pca_fit(intensity_mat, ntop = Inf, cent = TRUE, scal = FALSE, scores_metadata = cohort)

pca_scores <- pca_model_info$scores %>%
  dplyr::select(
    "SampleID" = "SampleID",
    "PC1" = ".fittedPC1",
    "PC2" = ".fittedPC2"
  )
write.csv(pca_scores, file = file.path(pcadir, "pca_scores12.csv"), row.names = FALSE)


saveRDS(pca_model_info, file = file.path(pcadir, "pca_model_info-all-vars.rds"))

message("PCA plots")

for (pca_colour_by in settings[["pca_colour_by"]]) {
  if (pca_colour_by == "Condition") {
    color_values <- colours_condition
  } else {
    color_values <- NULL
  }
  pc_12_plt <- pc_scores_loadings_plot(
    pca_model_info,
    pcx = 1, pcy = 2,
    color = pca_colour_by,
    color_values = color_values,
    loading_label = "Protein"
  )
  saveRDS(pc_12_plt, file = file.path(pcadir, sprintf("pca_12_%s-all-vars.rds", pca_colour_by)))
  ggsave(
    filename = file.path(pcadir, sprintf("pca_12_%s-all-vars.pdf", pca_colour_by)),
    plot = pc_12_plt,
    dpi = 300, width = unit(10, "in")
  ) # , height = unit(5, "in"))

  pc_34_plt <- pc_scores_loadings_plot(
    pca_model_info,
    pcx = 3, pcy = 4,
    color = pca_colour_by, color_values = color_values,
    loading_label = "Protein"
  )

  saveRDS(pc_34_plt, file = file.path(pcadir, sprintf("pca_34_%s-all-vars.rds", pca_colour_by)))
  ggsave(
    filename = file.path(pcadir, sprintf("pca_34_%s-all-vars.pdf", pca_colour_by)),
    plot = pc_34_plt,
    dpi = 300, width = unit(10, "in") # , height = unit(5, "in")
  )
}