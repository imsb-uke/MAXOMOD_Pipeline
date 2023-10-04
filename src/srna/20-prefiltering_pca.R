#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/pca_fit.R"))
source(here::here("src/shared/R-helpers/pc_scores_loadings_plot.R"))

suppressPackageStartupMessages({
  library("cowplot")
  library("dplyr")
  library("DESeq2")
  library("forcats")
  library("ggplot2")
  library("matrixStats")
  library("pheatmap")
  library("purrr")
  library("RColorBrewer")
  library("readr")
  library("tibble")
  library("tidyr")
})

settings <- get_settings(omic="srna", stage="prefiltering_pca", dataset_name=NULL) # dataset_name="c9"

cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = cols()
) %>% filter(! SampleID %in% settings[["exclude_samples"]])

disease_grp <- settings[["disease"]]
colours_sex <- c("female" = "#df65b0", "male" = "#2c7fb8")
colours_condition <- c("ctrl" = "#fe9929", "mut" = "#31a354")
names(colours_condition) <- c("ctrl", disease_grp)


dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

write.csv(cohort, file.path(settings[["output_directory"]], "cohort_excluded_outliers.csv"), row.names = FALSE)


# Read counts:
hairpin_filelist <- as.character(
  fs::dir_ls(path = file.path(settings[["input_directory"]], "results/bowtie/miRBase_hairpin"),
             glob = '*.hairpin.stats'
  )
)
names(hairpin_filelist) <- gsub(pattern = ".hairpin.stats$", "", basename(hairpin_filelist))

mature_filelist <- as.character(
  fs::dir_ls(path = file.path(settings[["input_directory"]], "results/bowtie/miRBase_mature"),
             glob = '*.mature.stats'
  )
)
names(mature_filelist) <- gsub(pattern = ".mature.stats$", "", basename(mature_filelist))

hairpin_counts <- map_dfr(hairpin_filelist, function(filename, sampleID) {
  readr::read_tsv(filename, col_names = c("mirnaID", "_", "counts", "_2"), col_types = readr::cols()) %>%
    dplyr::select(mirnaID, counts)
}, .id = "SampleID") %>%
  pivot_wider(names_from = "SampleID", values_from = "counts") %>%
  dplyr::select(mirnaID, one_of(cohort$SampleID))

mature_counts <- map_dfr(mature_filelist, function(filename, sampleID) {
  readr::read_tsv(filename, col_names = c("mirnaID", "_", "counts", "_2"), col_types = readr::cols()) %>%
    dplyr::select(mirnaID, counts)
}, .id = "SampleID") %>%
  pivot_wider(names_from = "SampleID", values_from = "counts") %>%
  dplyr::select(mirnaID, one_of(cohort$SampleID))



write.csv(x = mature_counts, file = file.path(settings[["output_directory"]], "mature_counts.csv"), row.names = FALSE)
saveRDS(mature_counts, file = file.path(settings[["output_directory"]], "mature_counts.rds"))

write.csv(x = hairpin_counts, file = file.path(settings[["output_directory"]], "hairpin_counts.csv"), row.names = FALSE)
saveRDS(hairpin_counts, file = file.path(settings[["output_directory"]], "hairpin_counts.rds"))



mature_counts_mat <- mature_counts %>%
  tibble::column_to_rownames("mirnaID")

hairpin_counts_mat <- hairpin_counts %>%
  tibble::column_to_rownames("mirnaID")

saveRDS(mature_counts_mat, file = file.path(settings[["output_directory"]], "mature_counts_mat.rds"))
saveRDS(hairpin_counts_mat, file = file.path(settings[["output_directory"]], "hairpin_counts_mat.rds"))


#hist(rowMeans(mature_counts_mat), breaks = c(0, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 300000),
#     xlim = c(0, 400), ylim = c(0, 0.01))

# Prepare count matrix

message(sprintf("Count matrix (%d mature and %d samples)", nrow(mature_counts_mat), ncol(mature_counts_mat)))
message(sprintf("Count matrix (%d hairpin and %d samples)", nrow(hairpin_counts_mat), ncol(hairpin_counts_mat)))



# Keep mirna that have more than ten counts in at least three samples of one group (Condition & Sex):

cohort_by_factors <- split(
  cohort,
  as.list(
    cohort[
      settings[["filter_low_counts"]][["stratify_by"]]
    ]
  )
)

mirna_mature_to_keep_list <- purrr::map(
  cohort_by_factors,
  function(cohort_i, matr, min_counts, at_least_samples, at_least_samples_mode) {
    if (at_least_samples_mode == "fraction") {
      min_samples <- ceiling(at_least_samples*nrow(cohort_i))
    } else if (at_least_samples_mode == "num_samples") {
      min_samples <- at_least_samples
    } else {
      stop("at_least_samples_mode must be either fraction or num_samples")
    }
    rnaseq_int_mat_i <- matr[,cohort_i$SampleID]
    ds <- DESeqDataSetFromMatrix(countData = rnaseq_int_mat_i,
                                 colData = cohort_i,
                                 design = ~ 1)
    keep <- rowSums(counts(ds) >= min_counts) >= min_samples
    names(keep[keep == TRUE])
  },
  matr = mature_counts_mat,
  min_counts = settings[["filter_low_counts"]][["at_least_counts"]],
  at_least_samples = settings[["filter_low_counts"]][["at_least_samples"]],
  at_least_samples_mode = settings[["filter_low_counts"]][["at_least_samples_mode"]]
)

mirna_hairpin_to_keep_list <- purrr::map(
  cohort_by_factors,
  function(cohort_i, matr, min_counts, at_least_samples, at_least_samples_mode) {
    if (at_least_samples_mode == "fraction") {
      min_samples <- ceiling(at_least_samples*nrow(cohort_i))
    } else if (at_least_samples_mode == "num_samples") {
      min_samples <- at_least_samples
    } else {
      stop("at_least_samples_mode must be either fraction or num_samples")
    }

    rnaseq_int_mat_i <- matr[,cohort_i$SampleID]
    ds <- DESeqDataSetFromMatrix(countData = rnaseq_int_mat_i,
                                 colData = cohort_i,
                                 design = ~ 1)
    keep <- rowSums(counts(ds) >= min_counts) >= min_samples
    names(keep[keep == TRUE])
  },
  matr = hairpin_counts_mat,
  min_counts = settings[["filter_low_counts"]][["at_least_counts"]],
  at_least_samples = settings[["filter_low_counts"]][["at_least_samples"]],
  at_least_samples_mode = settings[["filter_low_counts"]][["at_least_samples_mode"]]
)


write(jsonlite::toJSON(mirna_mature_to_keep_list, pretty = TRUE), file.path(settings[["output_directory"]], "to_keep_stratified_mirna_mature.json"))

pdf(file=file.path(settings[["output_directory"]], "to_keep_upset_mirna_mature.pdf"),
    width = 7, height = 4)
print(UpSetR::upset(UpSetR::fromList(mirna_mature_to_keep_list),  order.by = "freq"), newpage = FALSE)
dev.off()

mirna_mature_to_keep <- unique(unlist(mirna_mature_to_keep_list))

write(jsonlite::toJSON(setNames(list(mirna_mature_to_keep), "to_keep_mirna_mature"), pretty = TRUE),
      file.path(settings[["output_directory"]], "to_keep_mirna_mature.json"))
saveRDS(mirna_mature_to_keep, file = file.path(settings[["output_directory"]], "to_keep_mirna_mature.rds"))

print(sprintf("%d mature mirna kept from %d mature mirna", length(mirna_mature_to_keep), nrow(mature_counts_mat)))



write(jsonlite::toJSON(mirna_hairpin_to_keep_list, pretty = TRUE), file.path(settings[["output_directory"]], "to_keep_stratified_mirna_hairpin.json"))

pdf(file=file.path(settings[["output_directory"]], "to_keep_upset_mirna_hairpin.pdf"),
    width = 7, height = 4)
print(UpSetR::upset(UpSetR::fromList(mirna_hairpin_to_keep_list),  order.by = "freq"), newpage = FALSE)
dev.off()

mirna_hairpin_to_keep <- unique(unlist(mirna_hairpin_to_keep_list))

write(jsonlite::toJSON(setNames(list(mirna_hairpin_to_keep), "to_keep_mirna_hairpin"), pretty = TRUE),
      file.path(settings[["output_directory"]], "to_keep_mirna_hairpin.json"))
saveRDS(mirna_hairpin_to_keep, file = file.path(settings[["output_directory"]], "to_keep_mirna_hairpin.rds"))

message(sprintf("%d hairpin mirna kept from %d hairpin mirna", length(mirna_hairpin_to_keep), nrow(hairpin_counts_mat)))





mature_dds <- DESeqDataSetFromMatrix(
  countData = mature_counts_mat[mirna_mature_to_keep,],
  colData = cohort,
  design = ~1
)

hairpin_dds <- DESeqDataSetFromMatrix(
  countData = hairpin_counts_mat[mirna_hairpin_to_keep,],
  colData = cohort,
  design = ~1
)


# Normalization

x <- preprocessCore::normalize.quantiles(x = assay(mature_dds, "counts"), copy = TRUE)
dimnames(x) <- dimnames(assay(mature_dds, "counts"))
assay(mature_dds, "quantile_norm_counts_log2") <- log2(x)
rm(x)
saveRDS(mature_dds, file = file.path(settings[["output_directory"]], "counts_norm_quantile_mature.rds"))
write.csv(
  assay(mature_dds, "quantile_norm_counts_log2"),
  file = file.path(settings[["output_directory"]], "counts_norm_quantile_mature.csv"),
  row.names=TRUE
)


x <- preprocessCore::normalize.quantiles(x = assay(hairpin_dds, "counts"), copy = TRUE)
dimnames(x) <- dimnames(assay(hairpin_dds, "counts"))
assay(hairpin_dds, "quantile_norm_counts_log2") <- log2(x)
rm(x)
saveRDS(hairpin_dds, file = file.path(settings[["output_directory"]], "counts_norm_quantile_hairpin.rds"))
write.csv(
  assay(hairpin_dds, "quantile_norm_counts_log2"), 
  file = file.path(settings[["output_directory"]], "counts_norm_quantile_hairpin.csv"), 
  row.names=TRUE
)



# PCA:

mature_pca_model_info <- pca_fit(t(assay(mature_dds, "quantile_norm_counts_log2")),
                                 ntop = 500, cent = TRUE, scal = TRUE, scores_metadata = cohort)

mature_pca_scores <- mature_pca_model_info$scores %>%
  dplyr::select(
    "SampleID" = "SampleID",
    "PC1" = ".fittedPC1",
    "PC2" = ".fittedPC2"
 )
write.csv(mature_pca_scores, file = file.path(settings[["output_directory"]], "mature_pca_scores12.csv"), row.names=FALSE)

saveRDS(mature_pca_model_info, file = file.path(settings[["output_directory"]], "mature_pca_model_info.rds"))

hairpin_pca_model_info <- pca_fit(t(assay(hairpin_dds, "quantile_norm_counts_log2")),
                                  ntop = 360, cent = TRUE, scal = TRUE, scores_metadata = cohort)

hairpin_pca_scores <- hairpin_pca_model_info$scores %>%
  dplyr::select(
    "SampleID" = "SampleID",
    "PC1" = ".fittedPC1",
    "PC2" = ".fittedPC2"
 )
write.csv(hairpin_pca_scores, file = file.path(settings[["output_directory"]], "hairpin_pca_scores12.csv"), row.names=FALSE)


saveRDS(hairpin_pca_model_info, file = file.path(settings[["output_directory"]], "hairpin_pca_model_info.rds"))

plot(
  x = seq_along(mature_pca_model_info$var_percent),
  y = cumsum(mature_pca_model_info$var_percent),
  ylab = "Cummulated variance (%)",
  xlab = "Num of Principal Components"
)

plot(
  x = seq_along(hairpin_pca_model_info$var_percent),
  y = cumsum(hairpin_pca_model_info$var_percent),
  ylab = "Cummulated variance (%)",
  xlab = "Num of Principal Components"
)


mature_pc_12_plt <- pc_scores_loadings_plot(mature_pca_model_info, pcx = 1, pcy = 2,
                                            color = "Condition", color_values = colours_condition,
                                            loading_label = "mirnaID")
saveRDS(mature_pc_12_plt, file = file.path(settings[["output_directory"]], "mature_pc_12.rds"))
ggsave(filename = file.path(settings[["output_directory"]], "mature_pc_12.pdf"),
       plot = mature_pc_12_plt,
       dpi=300, width = unit(6, "in"), height = unit(5, "in"))


mature_pc_34_plt <- pc_scores_loadings_plot(mature_pca_model_info, pcx = 3, pcy = 4,
                                            color = "Condition", color_values = colours_condition,
                                            loading_label = "mirnaID")

saveRDS(mature_pc_34_plt, file = file.path(settings[["output_directory"]], "mature_pc_34.rds"))
ggsave(filename = file.path(settings[["output_directory"]], "mature_pc_34.pdf"),
       plot = mature_pc_34_plt,
       dpi=300, width = unit(6, "in"), height = unit(5, "in"))


hairpin_pc_12_plt <- pc_scores_loadings_plot(hairpin_pca_model_info, pcx = 1, pcy = 2,
                                             color = "Condition", color_values = colours_condition,
                                             loading_label = "mirnaID")
saveRDS(hairpin_pc_12_plt, file = file.path(settings[["output_directory"]], "hairpin_pc_12.rds"))
ggsave(filename = file.path(settings[["output_directory"]], "hairpin_pc_12.pdf"),
       plot = hairpin_pc_12_plt,
       dpi=300, width = unit(6, "in"), height = unit(5, "in"))


hairpin_pc_34_plt <- pc_scores_loadings_plot(hairpin_pca_model_info, pcx = 3, pcy = 4,
                                             color = "Condition", color_values = colours_condition,
                                             loading_label = "mirnaID")

saveRDS(hairpin_pc_34_plt, file = file.path(settings[["output_directory"]], "hairpin_pc_34.rds"))
ggsave(filename = file.path(settings[["output_directory"]], "hairpin_pc_34.pdf"),
       plot = hairpin_pc_34_plt,
       dpi=300, width = unit(6, "in"), height = unit(5, "in"))


if ("Center" %in% colnames(cohort)) {
  mature_pc_12_plt <- pc_scores_loadings_plot(mature_pca_model_info, pcx = 1, pcy = 2,
                                              color = "Center",
                                              loading_label = "mirnaID")
  saveRDS(mature_pc_12_plt, file = file.path(settings[["output_directory"]], "mature_pc_12-center.rds"))
  ggsave(filename = file.path(settings[["output_directory"]], "mature_pc_12-center.pdf"),
         plot = mature_pc_12_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))


  mature_pc_34_plt <- pc_scores_loadings_plot(mature_pca_model_info, pcx = 3, pcy = 4,
                                              color = "Center",
                                              loading_label = "mirnaID")

  saveRDS(mature_pc_34_plt, file = file.path(settings[["output_directory"]], "mature_pc_34-center.rds"))
  ggsave(filename = file.path(settings[["output_directory"]], "mature_pc_34-center.pdf"),
         plot = mature_pc_34_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))


  hairpin_pc_12_plt <- pc_scores_loadings_plot(hairpin_pca_model_info, pcx = 1, pcy = 2,
                                               color = "Center",
                                               loading_label = "mirnaID")
  saveRDS(hairpin_pc_12_plt, file = file.path(settings[["output_directory"]], "hairpin_pc_12-center.rds"))
  ggsave(filename = file.path(settings[["output_directory"]], "hairpin_pc_12-center.pdf"),
         plot = hairpin_pc_12_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))


  hairpin_pc_34_plt <- pc_scores_loadings_plot(hairpin_pca_model_info, pcx = 3, pcy = 4,
                                               color = "Center",
                                               loading_label = "mirnaID")

  saveRDS(hairpin_pc_34_plt, file = file.path(settings[["output_directory"]], "hairpin_pc_34-center.rds"))
  ggsave(filename = file.path(settings[["output_directory"]], "hairpin_pc_34-center.pdf"),
         plot = hairpin_pc_34_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))


}


#### PCA with all vars

mature_pca_model_info <- pca_fit(t(assay(mature_dds, "quantile_norm_counts_log2")),
                                 ntop = Inf, cent = TRUE, scal = TRUE, scores_metadata = cohort)
saveRDS(mature_pca_model_info, file = file.path(settings[["output_directory"]], "mature_pca_model_info-all-vars.rds"))

hairpin_pca_model_info <- pca_fit(t(assay(hairpin_dds, "quantile_norm_counts_log2")),
                                  ntop = Inf, cent = TRUE, scal = TRUE, scores_metadata = cohort)
saveRDS(hairpin_pca_model_info, file = file.path(settings[["output_directory"]], "hairpin_pca_model_info-all-vars.rds"))

plot(
  x = seq_along(mature_pca_model_info$var_percent),
  y = cumsum(mature_pca_model_info$var_percent),
  ylab = "Cummulated variance (%)",
  xlab = "Num of Principal Components"
)

plot(
  x = seq_along(hairpin_pca_model_info$var_percent),
  y = cumsum(hairpin_pca_model_info$var_percent),
  ylab = "Cummulated variance (%)",
  xlab = "Num of Principal Components"
)


mature_pc_12_plt <- pc_scores_loadings_plot(mature_pca_model_info, pcx = 1, pcy = 2,
                                            color = "Condition", color_values = colours_condition,
                                            loading_label = "mirnaID")
saveRDS(mature_pc_12_plt, file = file.path(settings[["output_directory"]], "mature_pc_12-all-vars.rds"))
ggsave(filename = file.path(settings[["output_directory"]], "mature_pc_12-all-vars.pdf"),
       plot = mature_pc_12_plt,
       dpi=300, width = unit(6, "in"), height = unit(5, "in"))


mature_pc_34_plt <- pc_scores_loadings_plot(mature_pca_model_info, pcx = 3, pcy = 4,
                                            color = "Condition", color_values = colours_condition,
                                            loading_label = "mirnaID")

saveRDS(mature_pc_34_plt, file = file.path(settings[["output_directory"]], "mature_pc_34-all-vars.rds"))
ggsave(filename = file.path(settings[["output_directory"]], "mature_pc_34-all-vars.pdf"),
       plot = mature_pc_34_plt,
       dpi=300, width = unit(6, "in"), height = unit(5, "in"))


hairpin_pc_12_plt <- pc_scores_loadings_plot(hairpin_pca_model_info, pcx = 1, pcy = 2,
                                             color = "Condition", color_values = colours_condition,
                                             loading_label = "mirnaID")
saveRDS(hairpin_pc_12_plt, file = file.path(settings[["output_directory"]], "hairpin_pc_12-all-vars.rds"))
ggsave(filename = file.path(settings[["output_directory"]], "hairpin_pc_12-all-vars.pdf"),
       plot = hairpin_pc_12_plt,
       dpi=300, width = unit(6, "in"), height = unit(5, "in"))


hairpin_pc_34_plt <- pc_scores_loadings_plot(hairpin_pca_model_info, pcx = 3, pcy = 4,
                                             color = "Condition", color_values = colours_condition,
                                             loading_label = "mirnaID")

saveRDS(hairpin_pc_34_plt, file = file.path(settings[["output_directory"]], "hairpin_pc_34-all-vars.rds"))
ggsave(filename = file.path(settings[["output_directory"]], "hairpin_pc_34-all-vars.pdf"),
       plot = hairpin_pc_34_plt,
       dpi=300, width = unit(6, "in"), height = unit(5, "in"))


if ("Center" %in% colnames(cohort)) {
  mature_pc_12_plt <- pc_scores_loadings_plot(mature_pca_model_info, pcx = 1, pcy = 2,
                                              color = "Center",
                                              loading_label = "mirnaID")
  saveRDS(mature_pc_12_plt, file = file.path(settings[["output_directory"]], "mature_pc_12-center-all-vars.rds"))
  ggsave(filename = file.path(settings[["output_directory"]], "mature_pc_12-center-all-vars.pdf"),
         plot = mature_pc_12_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))


  mature_pc_34_plt <- pc_scores_loadings_plot(mature_pca_model_info, pcx = 3, pcy = 4,
                                              color = "Center",
                                              loading_label = "mirnaID")

  saveRDS(mature_pc_34_plt, file = file.path(settings[["output_directory"]], "mature_pc_34-center-all-vars.rds"))
  ggsave(filename = file.path(settings[["output_directory"]], "mature_pc_34-center-all-vars.pdf"),
         plot = mature_pc_34_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))


  hairpin_pc_12_plt <- pc_scores_loadings_plot(hairpin_pca_model_info, pcx = 1, pcy = 2,
                                               color = "Center",
                                               loading_label = "mirnaID")
  saveRDS(hairpin_pc_12_plt, file = file.path(settings[["output_directory"]], "hairpin_pc_12-center-all-vars.rds"))
  ggsave(filename = file.path(settings[["output_directory"]], "hairpin_pc_12-center-all-vars.pdf"),
         plot = hairpin_pc_12_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))


  hairpin_pc_34_plt <- pc_scores_loadings_plot(hairpin_pca_model_info, pcx = 3, pcy = 4,
                                               color = "Center",
                                               loading_label = "mirnaID")

  saveRDS(hairpin_pc_34_plt, file = file.path(settings[["output_directory"]], "hairpin_pc_34-center-all-vars.rds"))
  ggsave(filename = file.path(settings[["output_directory"]], "hairpin_pc_34-center-all-vars.pdf"),
         plot = hairpin_pc_34_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))


}


# Heatmap

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


sampleDists <- dist(t(assay(mature_dds, "quantile_norm_counts_log2")))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix),
                                    mature_dds$Condition,
                                    substr(mature_dds$Sex,1,1),
                                    sep = "_")
colnames(sampleDistMatrix) <- NULL
plt <- pheatmap::pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
)

save_pheatmap_pdf(plt, file.path(settings[["output_directory"]], "mature_sample_heatmap.pdf"))


sampleDists <- dist(t(assay(hairpin_dds, "quantile_norm_counts_log2")))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix),
                                    hairpin_dds$Condition,
                                    substr(hairpin_dds$Sex,1,1),
                                    sep = "_")
colnames(sampleDistMatrix) <- NULL
plt <- pheatmap::pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
)

save_pheatmap_pdf(plt, file.path(settings[["output_directory"]], "hairpin_sample_heatmap.pdf"))


message("DONE")
