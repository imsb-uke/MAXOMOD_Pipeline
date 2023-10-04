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


settings <- get_settings(omic = "proteomics", stage = "prefiltering_pca", dataset_name = NULL) #dataset_name = "sod1"

cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = NULL
) %>% dplyr::filter(!SampleID %in% settings[["exclude_samples"]])


disease_grp <- setdiff(cohort$Condition, "ctrl")
colours_sex <- c("female" = "#df65b0", "male" = "#2c7fb8")
colours_condition <- c("ctrl" = "#fe9929", "als" = "#31a354")
names(colours_condition) <- c("ctrl", disease_grp)


dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)


protein_intensity <- readRDS(file.path(settings[["input_directory"]], "intensity.rds"))
protein_annotations <- readRDS(file.path(settings[["input_directory"]], "protein_annotations.rds"))

intensity_mat <- protein_intensity %>%
  tibble::column_to_rownames("UniProtName") %>%
  as.matrix() %>%
  t()

if (!all(rowSums(is.na(intensity_mat)) == 0)) {
  stop("Missing values detected. Unexpected.")
}


gplt <- tibble::enframe(sort(rowSums(intensity_mat == 0), decreasing = TRUE),
                        name = "SampleID",
                        value = "CountofZeroValues") %>%
  left_join(cohort, by = "SampleID") %>%
  ggplot() +
  geom_point(aes(x = forcats::fct_rev(forcats::fct_reorder(SampleID, -CountofZeroValues)),
                 y = CountofZeroValues,
                 color = Sex,
                 shape = Condition),
             size = 3) +
  labs(
    x = "Sample ID", y = "Num. proteins with zero intensity",
    title = "Distribution of zero intensity values across samples",
    subtitle = "GOOD: No pattern visible"
  ) +
  scale_color_manual(values = colours_sex) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()

ggsave(filename = file.path(settings[["output_directory"]], "zero_intensity_across_samples.pdf"), plot = gplt, height = unit(4, "in"), width = unit(6, "in"))
saveRDS(gplt, file.path(settings[["output_directory"]], "zero_intensity_across_samples.rds"))

intensity_mat <- intensity_mat[cohort$SampleID,]
protein_intensity <- protein_intensity %>% dplyr::select(UniProtName, one_of(cohort$SampleID))


# Pre-filtering: Keep samples with < n missing values in proteins
message("Pre-filtering: Remove samples with a lot of missing protein values")


if (!is.null(settings[["min_proteins_per_sample"]])) {
  samples_to_keep_nafilt <- rowSums(intensity_mat > 0) >= settings[["min_proteins_per_sample"]]
  samples_to_keep_nafilt <- names(samples_to_keep_nafilt)[samples_to_keep_nafilt]
  removed_samples <- setdiff(rownames(intensity_mat), samples_to_keep_nafilt)
  write(jsonlite::toJSON(removed_samples, pretty = TRUE),
        file.path(settings[["output_directory"]], "removed_samples_with_high_missing_values.json"))
  intensity_mat <- intensity_mat[samples_to_keep_nafilt,]
  protein_intensity <- protein_intensity %>% dplyr::select(UniProtName, one_of(samples_to_keep_nafilt))
  cohort <- cohort %>% dplyr::filter(SampleID %in% samples_to_keep_nafilt)
}


gplt <- tibble::enframe(sort(rowSums(intensity_mat == 0), decreasing = TRUE),
                        name = "SampleID",
                        value = "CountofZeroValues") %>%
  left_join(cohort, by = "SampleID") %>%
  ggplot() +
  geom_point(aes(x = forcats::fct_rev(forcats::fct_reorder(SampleID, -CountofZeroValues)),
                 y = CountofZeroValues,
                 color = Sex,
                 shape = Condition),
             size = 3) +
  labs(
    x = "Sample ID", y = "Num. proteins with zero intensity",
    title = "Distribution of zero intensity values across samples",
    subtitle = "GOOD: No pattern visible"
  ) +
  scale_color_manual(values = colours_sex) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()

ggsave(filename = file.path(settings[["output_directory"]], "zero_intensity_across_samples_filtered.pdf"), plot = gplt, height = unit(4, "in"), width = unit(6, "in"))
saveRDS(gplt, file.path(settings[["output_directory"]], "zero_intensity_across_samples_filtered.rds"))

message("Showing samples with lots of zero values")

x <- colSums(intensity_mat == 0)
proteins_with_some_zero_intensity <- tibble::enframe(sort(x[x > 0], decreasing = TRUE), name = "UniProtName", value = "CountofZeroValues")
proteins_with_some_zero_intensity

message(paste(c("A total of", length(x[x>0]), "out of", length(x), "proteins have zero intensity values"), sep = "", collapse =" "))


x <- as.data.frame(table(proteins_with_some_zero_intensity$CountofZeroValues), stringsAsFactors = FALSE)
x %>% dplyr::select("Number of proteins" = "Freq", "Samples with zero values" = "Var1") %>%
  dplyr::arrange(desc(as.numeric(`Samples with zero values`)))


# Pre-filtering: Remove proteins that have zero intensity in 3 or more samples of a given Sex&Condition
message("Pre-filtering proteins with low intensities")

if (!is.null(settings[["filter_low_counts"]][["stratify_by"]])) {
  cohort_by_factors <- split(
    cohort,
    as.list(
      cohort[
        settings[["filter_low_counts"]][["stratify_by"]]
      ]
    )
  )
} else {
  cohort_by_factors <- list(all=cohort)
}

protein_to_keep_list <- purrr::map(
  cohort_by_factors,
  function(cohort_i, matr, min_counts, at_least_samples, at_least_samples_mode) {
    if (at_least_samples_mode == "fraction") {
      min_samples <- ceiling(at_least_samples*nrow(cohort_i))
    } else if (at_least_samples_mode == "num_samples") {
      min_samples <- at_least_samples
    } else {
      stop("at_least_samples_mode must be either fraction or num_samples")
    }

    prot_mat_i <- matr[cohort_i$SampleID, ]

    keep <- colSums(prot_mat_i >= min_counts) >= min_samples
    names(keep[keep == TRUE])
  },
  matr = intensity_mat,
  min_counts = settings[["filter_low_counts"]][["at_least_counts"]],
  at_least_samples = settings[["filter_low_counts"]][["at_least_samples"]],
  at_least_samples_mode = settings[["filter_low_counts"]][["at_least_samples_mode"]]
)


write(jsonlite::toJSON(protein_to_keep_list, pretty = TRUE), file.path(settings[["output_directory"]], "to_keep_stratified_proteomics.json"))

if (length(protein_to_keep_list) > 1) {
  UpSetR::upset(UpSetR::fromList(protein_to_keep_list),  order.by = "freq")
  
  pdf(file=file.path(settings[["output_directory"]], "to_keep_stratified_proteomics.pdf"),
      width = 7, height = 4)
  print(UpSetR::upset(UpSetR::fromList(protein_to_keep_list),  order.by = "freq"), newpage = FALSE)
  dev.off()
}

protein_to_keep <- unique(unlist(protein_to_keep_list))


write(jsonlite::toJSON(setNames(list(protein_to_keep), "protein_to_keep"), pretty = TRUE),
      file.path(settings[["output_directory"]], "protein_to_keep.json"))
saveRDS(protein_to_keep, file = file.path(settings[["output_directory"]], "protein_to_keep.json"))

message("Kept ", length(protein_to_keep), " proteins")

protein_intensity <- protein_intensity %>%
  dplyr::filter(UniProtName %in% protein_to_keep)

intensity_mat <- intensity_mat[,colnames(intensity_mat) %in% protein_to_keep]


write.csv(x = intensity_mat, file = file.path(settings[["output_directory"]], "intensity_mat_filtered.csv"), row.names = TRUE)
saveRDS(intensity_mat, file = file.path(settings[["output_directory"]], "intensity_mat_filtered.rds"))



# We can represent the intensities of all the samples for each of those proteins to inspect their distribution:
message("Some plots...")

x <- protein_intensity %>%
  pivot_longer(-UniProtName, names_to = "SampleID", values_to = "Intensity") %>%
  left_join(cohort, by = "SampleID")

ggplot(x) +
  geom_boxplot(aes(x = forcats::fct_reorder(SampleID, as.numeric(as.factor(Condition))), y = Intensity, fill = Sex, color = Condition)) +
  labs(x = "SampleID", title = "LFQ intensities") +
  scale_y_sqrt() +
  scale_fill_manual(values = colours_sex) +
  scale_colour_manual(values = colours_condition) +
  coord_flip()


if ("XPO1_MOUSE" %in% x$UniProtName) {
  ggplot(dplyr::filter(x, UniProtName == "XPO1_MOUSE")) +
    geom_point(aes(x = forcats::fct_reorder(SampleID, Condition), y = Intensity, shape = Sex, color = Condition), size = 3) +
    labs(x = "SampleID", title = "XPO1_MOUSE raw intensity") +
    scale_y_sqrt() +
    scale_colour_manual(values = colours_condition) +
    coord_flip()
}

if ("XPO1_HUMAN" %in% x$UniProtName) {
  ggplot(dplyr::filter(x, UniProtName == "XPO1_HUMAN")) +
    geom_point(aes(x = forcats::fct_reorder(SampleID, Condition), y = Intensity, shape = Sex, color = Condition), size = 3) +
    labs(x = "SampleID", title = "XPO1_HUMAN raw intensity") +
    scale_y_sqrt() +
    scale_colour_manual(values = colours_condition) +
    coord_flip()
}




# Imputation:
message("Imputation")

if (settings[["zero_value_imputation"]][["method"]] == "min_observed") {
  min_intensity <- min(intensity_mat[intensity_mat > 0])
  intensity_mat[intensity_mat == 0] <- min_intensity
  protein_intensity <- protein_intensity %>%
    mutate(across(-UniProtName, function(x) { x[x==0] <- min_intensity; x }))
  rm(min_intensity)
} else if (settings[["zero_value_imputation"]][["method"]] == "min_observed_in_protein") {
  min_intensities <- protein_intensity %>%
    tidyr::pivot_longer(cols=-UniProtName, names_to = "SampleID", values_to = "intensity") %>%
    dplyr::group_by(UniProtName) %>%
    dplyr::filter(intensity > 0) %>%
    dplyr::summarise(min_intensity = min(intensity)) %>%
    dplyr::ungroup() %>%
    tibble::deframe()
  stopifnot(all(c("UniProtName", rownames(intensity_mat)) == colnames(protein_intensity)))
  for (protein in colnames(intensity_mat)) {
    protein_values <- intensity_mat[,protein]
    protein_values[protein_values==0] <- min_intensities[protein]
    intensity_mat[,protein] <- protein_values
    protein_intensity[protein_intensity$UniProtName == protein, rownames(intensity_mat)] <- as.list(protein_values)
  }
  rm(protein_values)
} else if (settings[["zero_value_imputation"]][["method"]] == "missForest") {
  intensity_mat[intensity_mat==0] <- NA
  misforout <- missForest::missForest(intensity_mat)
  intensity_mat <- misforout$ximp
  protein_intensity <- t(intensity_mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("UniProtName") %>%
    dplyr::select(one_of(colnames(protein_intensity)))
} else if (settings[["zero_value_imputation"]][["method"]] == "none") {

} else {
  stop("Invalid imputation method")
}


write.csv(x = intensity_mat, file = file.path(settings[["output_directory"]], "intensity_mat_filtered_imputed.csv"), row.names = TRUE)
saveRDS(intensity_mat, file = file.path(settings[["output_directory"]], "intensity_mat_filtered_imputed.rds"))


# Log2 transform
message("Log2 transform")

protein_intensity_log <- protein_intensity %>%
  mutate(across(-UniProtName, log2))
protein_intensity_log[1:20, 1:3]

intensity_log_mat <- log2(intensity_mat)

write.csv(x = intensity_log_mat, file = file.path(settings[["output_directory"]], "intensity_mat_filtered_imputed_log2transf.csv"), row.names = TRUE)
saveRDS(intensity_log_mat, file = file.path(settings[["output_directory"]], "intensity_mat_filtered_imputed_log2transf.rds"))


# We can represent the intensities of all the samples for each of those proteins to inspect their distribution:

x <- protein_intensity_log %>%
  pivot_longer(-UniProtName, names_to = "SampleID", values_to = "Intensity") %>%
  left_join(cohort, by = "SampleID")

ggplot(x) +
  geom_boxplot(aes(x = forcats::fct_reorder(SampleID, as.numeric(as.factor(Condition))), y = Intensity, fill = Sex, color = Condition)) +
  labs(x = "SampleID", title = "raw intensities") +
  scale_fill_manual(values = colours_sex) +
  scale_colour_manual(values = colours_condition) +
  coord_flip()


if ("XPO1_MOUSE" %in% x$UniProtName) {
  ggplot(filter(x, UniProtName == "XPO1_MOUSE")) +
    geom_point(aes(x = forcats::fct_reorder(SampleID, Condition), y = Intensity, shape = Sex, color = Condition), size = 3) +
    labs(x = "SampleID", title = "XPO1_MOUSE log2transf intensity") +
    scale_colour_manual(values = colours_condition) +
    coord_flip()
}



# Normalization
message("Normalization")

# norm_pqn <- function(spectra) {
#   num_samples <- nrow(spectra)
#   if (num_samples < 10) {
#     warning("The Probabalistic Quotient Normalization requires several samples ",
#             "to compute the median spectra. Your number of samples is low")
#   }
#   # Normalize to the area
#   areas <- rowSums(spectra)
#   spectra2 <- spectra / areas
#   if (num_samples == 1) {
#     # We have warned, and here there is nothing to do anymore
#     warning("PQN is absurd with a single sample. We have normalized it to the area.")
#     return(list(spectra = spectra2,
#                 norm_factor = areas))
#   }
#   # Move spectra above zero:
#   if (any(spectra2 < 0)) {
#     spectra2 <- spectra2 - min(spectra2)
#   }
#   # Median of each protein: (We need multiple spectra in order to get a reliable median!)
#   m <- matrixStats::colMedians(spectra2)
#   # Divide at each ppm by its median:
#   f <- spectra2/m[col(spectra2)]
#   f <- matrixStats::rowMedians(f)
#   # Divide each spectra by its f value
#   spectra3 <- spectra2 / f
#   list(spectra = spectra3,
#        norm_factor = f*areas)
# }
#

#intensity_log_mat_norm2 <- t(vsn::justvsn(t(intensity_log_mat)))
#intensity_log_mat_norm_detail <- norm_pqn(intensity_log_mat)
#intensity_log_mat_norm <- intensity_log_mat_norm_detail$spectra

norm_method <- as.character(settings[["normalization"]][["method"]])
if (norm_method == "quantiles") {
  intensity_log_mat_norm <- t(preprocessCore::normalize.quantiles(t(intensity_log_mat),copy=TRUE))
  rownames(intensity_log_mat_norm) <- rownames(intensity_log_mat)
  colnames(intensity_log_mat_norm) <- colnames(intensity_log_mat)

  protein_intensity_log_norm <- tibble::as_tibble(t(intensity_log_mat_norm), rownames = "UniProtName")
} else if (norm_method == "none") {
  intensity_log_mat_norm <- intensity_log_mat
  protein_intensity_log_norm <- protein_intensity
} else {
  stop("Unknown normalization method")
}

write.csv(x = intensity_log_mat_norm, file = file.path(settings[["output_directory"]], "intensity_mat_filtered_imputed_log2transf_norm.csv"), row.names = TRUE)
saveRDS(intensity_log_mat_norm, file = file.path(settings[["output_directory"]], "intensity_mat_filtered_imputed_log2transf_norm.rds"))


x <- protein_intensity_log_norm %>%
  pivot_longer(-UniProtName, names_to = "SampleID", values_to = "Intensity") %>%
  left_join(cohort, by = "SampleID")


ggplot(x) +
  geom_boxplot(aes(x = forcats::fct_reorder(SampleID, as.numeric(as.factor(Condition))), y = Intensity, fill = Sex, color = Condition)) +
  labs(x = "SampleID", title = "raw intensities") +
  scale_fill_manual(values = colours_sex) +
  scale_colour_manual(values = colours_condition) +
  coord_flip()


if ("XPO1_MOUSE" %in% x$UniProtName) {
  ggplot(filter(x, UniProtName == "XPO1_MOUSE")) +
    geom_point(aes(x = forcats::fct_reorder(SampleID, Condition), y = Intensity, shape = Sex, color = Condition), size = 3) +
    labs(x = "SampleID", title = "XPO1_MOUSE normalized intensity") +
    scale_colour_manual(values = colours_condition) +
    coord_flip()
}


# PCA --------------
message("Creating PCA")
pca_model_info <- pca_fit(intensity_log_mat_norm, ntop = 500, cent = TRUE, scal = FALSE, scores_metadata = cohort)
saveRDS(pca_model_info, file = file.path(settings[["output_directory"]], "pca_model_info.rds"))

# plot(
#   x = seq_along(pca_model_info$var_percent),
#   y = cumsum(pca_model_info$var_percent),
#   ylab = "Cummulated variance (%)",
#   xlab = "Num of Principal Components"
# )

#table(cohort_for_pca$Condition, cohort_for_pca$Center)

message("PCA plots")


for (pca_colour_by in settings[["pca_colour_by"]]) {
  if (pca_colour_by == "Condition") {
    color_values <- colours_condition
  } else {
    color_values <- NULL
  }
  pc_12_plt <- pc_scores_loadings_plot(
    pca_model_info, pcx = 1, pcy = 2,
    color = pca_colour_by,
    color_values = color_values,
    loading_label = "Protein"
  )
  saveRDS(pc_12_plt, file = file.path(settings[["output_directory"]], sprintf("pca_12_%s.rds", pca_colour_by)))
  ggsave(filename = file.path(settings[["output_directory"]], sprintf("pca_12_%s.pdf", pca_colour_by)),
         plot = pc_12_plt,
         dpi=300, width = unit(10, "in"))#, height = unit(5, "in"))

  pc_34_plt <- pc_scores_loadings_plot(
    pca_model_info, pcx = 3, pcy = 4,
    color = pca_colour_by, color_values = color_values,
    loading_label = "Protein"
  )

  saveRDS(pc_34_plt, file = file.path(settings[["output_directory"]], sprintf("pca_34_%s.rds", pca_colour_by)))
  ggsave(
    filename = file.path(settings[["output_directory"]], sprintf("pca_34_%s.pdf", pca_colour_by)),
    plot = pc_34_plt,
    dpi=300, width = unit(10, "in")#, height = unit(5, "in")
  )

}



# PCA - all vars --------------
message("Creating PCA")
pca_model_info <- pca_fit(intensity_log_mat_norm, ntop = Inf, cent = TRUE, scal = FALSE, scores_metadata = cohort)

pca_scores <- pca_model_info$scores %>%
  dplyr::select(
    "SampleID" = "SampleID",
    "PC1" = ".fittedPC1",
    "PC2" = ".fittedPC2"
 )
write.csv(pca_scores, file = file.path(settings[["output_directory"]], "pca_scores12.csv"), row.names=FALSE)


saveRDS(pca_model_info, file = file.path(settings[["output_directory"]], "pca_model_info-all-vars.rds"))

# plot(
#   x = seq_along(pca_model_info$var_percent),
#   y = cumsum(pca_model_info$var_percent),
#   ylab = "Cummulated variance (%)",
#   xlab = "Num of Principal Components"
# )

#table(cohort_for_pca$Condition, cohort_for_pca$Center)

message("PCA plots")



for (pca_colour_by in settings[["pca_colour_by"]]) {
  if (pca_colour_by == "Condition") {
    color_values <- colours_condition
  } else {
    color_values <- NULL
  }
  pc_12_plt <- pc_scores_loadings_plot(
    pca_model_info, pcx = 1, pcy = 2,
    color = pca_colour_by,
    color_values = color_values,
    loading_label = "Protein"
  )
  saveRDS(pc_12_plt, file = file.path(settings[["output_directory"]], sprintf("pca_12_%s-all-vars.rds", pca_colour_by)))
  ggsave(filename = file.path(settings[["output_directory"]], sprintf("pca_12_%s-all-vars.pdf", pca_colour_by)),
         plot = pc_12_plt,
         dpi=300, width = unit(10, "in"))#, height = unit(5, "in"))

  pc_34_plt <- pc_scores_loadings_plot(
    pca_model_info, pcx = 3, pcy = 4,
    color = pca_colour_by, color_values = color_values,
    loading_label = "Protein"
  )

  saveRDS(pc_34_plt, file = file.path(settings[["output_directory"]], sprintf("pca_34_%s-all-vars.rds", pca_colour_by)))
  ggsave(
    filename = file.path(settings[["output_directory"]], sprintf("pca_34_%s-all-vars.pdf", pca_colour_by)),
    plot = pc_34_plt,
    dpi=300, width = unit(10, "in")#, height = unit(5, "in")
  )

}


# Write cohort filtered
write.csv(x = cohort, file = file.path(settings[["output_directory"]], "cohort_filtered.csv"), row.names = FALSE)

