suppressPackageStartupMessages({
  requireNamespace("jsonlite")
  requireNamespace("DESeq2")
  requireNamespace("tibble")
  requireNamespace("dplyr")
})

deseq2_analysis <- function(
  count_data,
  col_data,
  deseq2_settings,
  output_directory,
  feature_column = "feature",
  results_callback = NULL
) {
  dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
  count_data_orig <- count_data
  col_data_orig <- col_data
  for (deseq2_setting in deseq2_settings) {
    deseq2_result_outdir <- file.path(output_directory, deseq2_setting[["name"]])
    dir.create(deseq2_result_outdir, showWarnings = FALSE, recursive = TRUE)
    count_data <- count_data_orig
    col_data <- col_data_orig
    if ("autoscale" %in% names(deseq2_setting)) {
      message("Autoscaling: ", paste(deseq2_setting[["autoscale"]], collapse = ", "))
      autoscale_stats <- list()
      for (scale_var in deseq2_setting[["autoscale"]]) {
        autoscale_stats[[scale_var]] <- list(
          mean = mean(col_data[[scale_var]], na.rm = TRUE),
          sd = sd(col_data[[scale_var]], na.rm = TRUE)
        )
        col_data[[scale_var]] <- scale(col_data[[scale_var]], center = TRUE, scale = TRUE)
      }
      jsonlite::write_json(autoscale_stats, path = file.path(deseq2_result_outdir, "autoscale_stats.json"))
    }
    if ("filter" %in% names(deseq2_setting)) {
      message("Filtering samples")
      for (criteria in deseq2_setting[["filter"]]) {
        samples_to_keep <- col_data[[criteria[["field"]]]] == criteria[["value"]]
        count_data <- count_data[, samples_to_keep]
        col_data <- col_data[samples_to_keep,]
      }
    }
    design_formula <- as.formula(deseq2_setting[["formula"]])
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = count_data,
      colData = col_data,
      design = design_formula
    )
    dds <- DESeq2::DESeq(dds)
    saveRDS(dds, file.path(deseq2_result_outdir, "dds.rds"))
    for (result_contrast in deseq2_setting[["results"]]) {
      r <- DESeq2::results(
        dds,
        contrast = as.list(result_contrast[["contrast"]]),
        pAdjustMethod = result_contrast[["padjust_method"]]
      )
      print(r)
      r <- as.data.frame(r)
      r <- tibble::rownames_to_column(r, feature_column)
      r <- dplyr::arrange(r, padj)
      if (!is.null(results_callback)) {
        r <- results_callback(r)
      }
      write.csv(
        x = r,
        file = file.path(deseq2_result_outdir, sprintf("%s.csv", result_contrast[["name"]])),
        row.names = FALSE
      )
      saveRDS(object = r, file = file.path(deseq2_result_outdir, sprintf("%s.rds", result_contrast[["name"]])))
      ## Histogram of p-values
      png(filename = file.path(deseq2_result_outdir, sprintf("%s_hist_pval.png", result_contrast[["name"]])))
      hist(r$pvalue[r$baseMean > 1], xlab = "p-values (genes with baseMean > 1)", main = "Histogram")
      dev.off()
    }
  }
}
