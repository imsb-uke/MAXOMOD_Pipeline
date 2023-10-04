#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/cohort_to_factors.R"))
source(here::here("src/shared/R-helpers/get_settings.R"))
source(here::here("src/shared/R-helpers/pca_fit.R"))
source(here::here("src/shared/R-helpers/pc_scores_loadings_plot.R"))


suppressPackageStartupMessages({
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



settings <- get_settings(omic = "rnaseq", stage = "prefiltering_pca", dataset_name = NULL) # dataset_name="sod1"

cohort <- readr::read_csv(
  settings[["cohort"]],
  col_types = NULL
) %>% dplyr::filter(!SampleID %in% settings[["exclude_samples"]])

if ("Center" %in% colnames(cohort)) {
  disease_grp <- "als"
} else {
  disease_grp <- "mut"
}

colours_sex <- c("female" = "#df65b0", "male" = "#2c7fb8")
colours_condition <- c("ctrl" = "#fe9929", "mut" = "#31a354")
names(colours_condition) <- c("ctrl", disease_grp)

cohort_nf_map <- readr::read_csv(
  settings[["nfid_to_sample_id"]],
  col_types = cols(
    SampleID = col_character(),
    nf_id = col_character()
  )
) %>%
  dplyr::select(nf_id, SampleID) %>%
  tibble::deframe()

cohort_nf_map[cohort_nf_map %in% cohort$SampleID]


# Read counts:
xx <- readRDS(settings[["counts"]])
xx <- xx[,rownames(colData(xx)) %in% names(cohort_nf_map[cohort_nf_map %in% cohort$SampleID])]
colData(xx)["SampleID"] <- cohort_nf_map[rownames(colData(xx))]
colData(xx) <- S4Vectors::DataFrame(tibble::column_to_rownames(dplyr::left_join(as.data.frame(colData(xx)), cohort, by = "SampleID"), "names"))

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

write.csv(cohort, file.path(settings[["output_directory"]], "cohort_excluded_outliers.csv"), row.names = FALSE)


# ENSG to symbol map:
ensg_symbol <- dplyr::select(as.data.frame(rowData(xx)), gene_id, gene_name)
utils::write.csv(ensg_symbol, file = file.path(settings[["output_directory"]], "ENSG_to_symbol.csv"), row.names = FALSE)
saveRDS(ensg_symbol, file = file.path(settings[["output_directory"]], "ENSG_to_symbol.rds"))


# Prepare count matrix
rnaseq_intensities_mat <- round(as.matrix(assays(xx)[["counts"]]))
mode(rnaseq_intensities_mat) <- "integer"
colnames(rnaseq_intensities_mat) <- cohort_nf_map[colnames(rnaseq_intensities_mat)]

write.csv(rnaseq_intensities_mat, file = file.path(settings[["output_directory"]], "counts_raw_all_samples.csv"), row.names = TRUE)
saveRDS(rnaseq_intensities_mat, file = file.path(settings[["output_directory"]], "counts_raw_all_samples.rds"))

rnaseq_intensities_mat <- rnaseq_intensities_mat[, cohort$SampleID]

print(sprintf("Count Matrix (%d genes and %d samples)", nrow(rnaseq_intensities_mat), ncol(rnaseq_intensities_mat)))
print(rnaseq_intensities_mat[1:5,1:5])

write.csv(rnaseq_intensities_mat, file = file.path(settings[["output_directory"]], "counts_raw_no_outliers.csv"), row.names = TRUE)
saveRDS(rnaseq_intensities_mat, file = file.path(settings[["output_directory"]], "counts_raw_no_outliers.rds"))

#tmp1 <- rowMeans(rnaseq_intensities_mat)
#ggplot2::qplot(x = tmp1, geom = "histogram", bins=50) + scale_x_log10()

# Counts for Xist:

if ("Xist" %in% ensg_symbol$gene_name) {
  xist_symbol <- "Xist"
  xist_symbol_sym <- rlang::sym(xist_symbol)
  xist_ensembl_id <- ensg_symbol %>% dplyr::filter(gene_name == xist_symbol) %>% pull("gene_id")
} else if ("XIST" %in% ensg_symbol$gene_name) {
  xist_symbol <- "XIST"
  xist_symbol_sym <- rlang::sym(xist_symbol)
  xist_ensembl_id <- ensg_symbol %>% dplyr::filter(gene_name == xist_symbol) %>% pull("gene_id")
} else {
  xist_symbol <- NULL
  xist_symbol_sym <- NULL
  xist_ensebl_id <- NULL
}


if (!is.null(xist_symbol)) {
  to_plot <- rnaseq_intensities_mat[xist_ensembl_id,] %>%
    tibble::enframe(name = "SampleID", value = xist_symbol) %>%
    dplyr::left_join(cohort, by = "SampleID") %>%
    dplyr::arrange(desc(!!xist_symbol_sym))

  gplt <- ggplot(to_plot) +
    geom_col(aes(x = fct_reorder(SampleID, !!xist_symbol_sym), y = !!xist_symbol_sym, fill= Sex)) +
    scale_colour_manual(values = colours_condition) +
    scale_fill_manual(values = colours_sex) +
    labs(x = "SampleID", y = sprintf("%s (raw counts)", xist_symbol)) +
    coord_flip()

  saveRDS(gplt, file = file.path(settings[["output_directory"]], sprintf("%s_counts_raw.rds", xist_symbol)))
  ggsave(filename = file.path(settings[["output_directory"]], sprintf("%s_counts_raw.pdf", xist_symbol)), plot = gplt,
         dpi=300, width = unit(6, "in"), height = unit(4, "in"))

}

# Keep genes that have more than ten counts in at least three samples of one group (Condition & Sex):

cohort_by_factors <- split(
  cohort,
  as.list(
    cohort[
      settings[["filter_low_counts"]][["stratify_by"]]
      ]
  )
)

genes_to_keep_list <- purrr::map(
  cohort_by_factors,
  function(cohort_i,
           rnaseq_intensities_mat,
           min_counts,
           at_least_samples,
           at_least_samples_mode
           ) {
    if (at_least_samples_mode == "fraction") {
      min_samples <- ceiling(at_least_samples*nrow(cohort_i))
    } else if (at_least_samples_mode == "num_samples") {
      min_samples <- at_least_samples
    } else {
      stop("at_least_samples_mode must be either fraction or num_samples")
    }
    rnaseq_int_mat_i <- rnaseq_intensities_mat[,cohort_i$SampleID]
    ds <- DESeqDataSetFromMatrix(countData = rnaseq_int_mat_i,
                                 colData = cohort_i,
                                 design = ~ 1)
    keep <- rowSums(counts(ds) >= min_counts) >= min_samples
    names(keep[keep == TRUE])
  },
  rnaseq_intensities_mat = rnaseq_intensities_mat,
  min_counts = settings[["filter_low_counts"]][["at_least_counts"]],
  at_least_samples = settings[["filter_low_counts"]][["at_least_samples"]],
  at_least_samples_mode = settings[["filter_low_counts"]][["at_least_samples_mode"]]
)

write(jsonlite::toJSON(genes_to_keep_list, pretty = TRUE), file.path(settings[["output_directory"]], "genes_to_keep_stratified.json"))

pdf(file=file.path(settings[["output_directory"]], "genes_to_keep_upset.pdf"),
    width = 7, height = 4)
print(UpSetR::upset(UpSetR::fromList(genes_to_keep_list),  order.by = "freq"), newpage = FALSE)
dev.off()

genes_to_keep <- unique(unlist(genes_to_keep_list))

write(jsonlite::toJSON(setNames(list(genes_to_keep), "genes_to_keep"), pretty = TRUE),
      file.path(settings[["output_directory"]], "genes_to_keep.json"))
saveRDS(genes_to_keep, file = file.path(settings[["output_directory"]], "genes_to_keep.rds"))


print(sprintf("%d genes kept from %d genes", length(genes_to_keep),nrow(rnaseq_intensities_mat)))

write.csv(rnaseq_intensities_mat, file = file.path(settings[["output_directory"]], "counts_raw_no_outliers.csv"), row.names = TRUE)
saveRDS(rnaseq_intensities_mat, file = file.path(settings[["output_directory"]], "counts_raw_no_outliers.rds"))

rnaseq_intensities_mat_kept <- rnaseq_intensities_mat[genes_to_keep,]
write.csv(
  rnaseq_intensities_mat_kept,
  file = file.path(settings[["output_directory"]], "counts_raw_no_outliers_no_low_genes.csv"),
  row.names = TRUE
)
saveRDS(
  rnaseq_intensities_mat_kept,
  file = file.path(settings[["output_directory"]], "counts_raw_no_outliers_no_low_genes.rds")
)

dds <- DESeqDataSetFromMatrix(
  countData = rnaseq_intensities_mat_kept,
  colData = cohort,
  design = ~1
)

# Normalization

vsd <- DESeq2::vst(dds, blind = TRUE)

write.csv(assay(vsd), file = file.path(settings[["output_directory"]], "counts_norm_vst.csv"), row.names = TRUE)
saveRDS(vsd, file = file.path(settings[["output_directory"]], "counts_norm_vst.rds"))

# Xist expression (normalized)

if (!is.null(xist_symbol)) {
  to_plot <- left_join(
    tibble::enframe(assay(vsd)[xist_ensembl_id,], name = "SampleID", value = xist_symbol),
    cohort,
    by = "SampleID"
  )

  gplt <- ggplot(to_plot) +
    geom_col(aes(x = fct_reorder(SampleID, !!xist_symbol_sym), y = !!xist_symbol_sym, fill= Sex)) +
    scale_colour_manual(values = colours_condition) +
    scale_fill_manual(values = colours_sex) +
    labs(x = "SampleID", y = sprintf("%s (normalized counts -- variance stabilized, blind design)", xist_symbol)) +
    coord_flip()

  saveRDS(gplt, file = file.path(settings[["output_directory"]], sprintf("%s_counts_vst.rds", xist_symbol)))
  ggsave(filename = file.path(settings[["output_directory"]], sprintf("%s_counts_vst.pdf", xist_symbol)), plot = gplt,
         dpi=300, width = unit(6, "in"), height = unit(4, "in"))

}


# PCA:

pca_model_info <- pca_fit(t(assay(vsd)), ntop = 500, cent = TRUE, scal = FALSE, scores_metadata = cohort)

pca_model_info$loadings <- pca_model_info$loadings %>%
  dplyr::left_join(ensg_symbol, by = c("column" = "gene_id")) %>%
  dplyr::select(column, gene_name, everything())



saveRDS(pca_model_info, file = file.path(settings[["output_directory"]], "pca_model_info.rds"))

# plot(
#   x = seq_along(pca_model_info$var_percent),
#   y = cumsum(pca_model_info$var_percent),
#   ylab = "Cummulated variance (%)",
#   xlab = "Num of Principal Components"
# )

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
    loading_colname = "gene_name",
    loading_label = "Gene"
  )
  saveRDS(pc_12_plt, file = file.path(settings[["output_directory"]], sprintf("pca_12_%s.rds", pca_colour_by)))
  ggsave(filename = file.path(settings[["output_directory"]], sprintf("pca_12_%s.pdf", pca_colour_by)),
         plot = pc_12_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))

  pc_34_plt <- pc_scores_loadings_plot(
    pca_model_info, pcx = 3, pcy = 4,
    color = pca_colour_by, color_values = color_values,
    loading_colname = "gene_name",
    loading_label = "Gene"
  )

  saveRDS(pc_34_plt, file = file.path(settings[["output_directory"]], sprintf("pca_34_%s.rds", pca_colour_by)))
  ggsave(
    filename = file.path(settings[["output_directory"]], sprintf("pca_34_%s.pdf", pca_colour_by)),
    plot = pc_34_plt,
    dpi=300, width = unit(6, "in"), height = unit(5, "in")
  )

}


pca_model_info <- pca_fit(t(assay(vsd)), ntop = Inf, cent = TRUE, scal = FALSE, scores_metadata = cohort)

pca_model_info$loadings <- pca_model_info$loadings %>%
  dplyr::left_join(ensg_symbol, by = c("column" = "gene_id")) %>%
  dplyr::select(column, gene_name, everything())


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
    loading_colname = "gene_name",
    loading_label = "Gene"
  )
  saveRDS(pc_12_plt, file = file.path(settings[["output_directory"]], sprintf("pca_12_%s-all-vars.rds", pca_colour_by)))
  ggsave(filename = file.path(settings[["output_directory"]], sprintf("pca_12_%s-all-vars.pdf", pca_colour_by)),
         plot = pc_12_plt,
         dpi=300, width = unit(6, "in"), height = unit(5, "in"))

  pc_34_plt <- pc_scores_loadings_plot(
    pca_model_info, pcx = 3, pcy = 4,
    color = pca_colour_by, color_values = color_values,
    loading_colname = "gene_name",
    loading_label = "Gene"
  )

  saveRDS(pc_34_plt, file = file.path(settings[["output_directory"]], sprintf("pca_34_%s-all-vars.rds", pca_colour_by)))
  ggsave(
    filename = file.path(settings[["output_directory"]], sprintf("pca_34_%s-all-vars.pdf", pca_colour_by)),
    plot = pc_34_plt,
    dpi=300, width = unit(6, "in"), height = unit(5, "in")
  )

}

# Heatmap

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )

if ("Center" %in% colnames(cohort)) {
  rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix),
                                      vsd$Condition,
                                      substr(vsd$Sex,1,1),
                                      vsd$Center,
                                      sep = "_")
} else {
  rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix),
                                      vsd$Condition,
                                      substr(vsd$Sex,1,1),
                                      sep = "_")

}

colnames(sampleDistMatrix) <- NULL
plt <- pheatmap::pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(plt, file.path(settings[["output_directory"]], "sample_heatmap.pdf"))
message("DONE")
