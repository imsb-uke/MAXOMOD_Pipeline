
requireNamespace("jsonlite")
requireNamespace("tibble")
requireNamespace("dplyr")
requireNamespace("limma")

limma_analysis <- function(
  x,
  col_data,
  limma_settings,
  output_directory,
  results_callback = NULL
) {
  # x has samples in rows, proteins in columns
  dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
  x_orig <- x
  col_data_orig <- col_data
  de_results <- list()
  for (analysis in limma_settings) {
    analysis_name <- analysis[["name"]]
    limma_result_outdir <- file.path(output_directory, analysis_name)
    dir.create(limma_result_outdir, showWarnings = FALSE, recursive = TRUE)
    x <- x_orig
    col_data <- col_data_orig
    if ("autoscale" %in% names(analysis)) {
      message("Autoscaling: ", paste(analysis[["autoscale"]], collapse = ", "))
      autoscale_stats <- list()
      for (scale_var in analysis[["autoscale"]]) {
        autoscale_stats[[scale_var]] <- list(
          mean= mean(col_data[[scale_var]], na.rm = TRUE),
          sd = sd(col_data[[scale_var]], na.rm = TRUE)
        )
        col_data[[scale_var]] <- scale(col_data[[scale_var]], center = TRUE, scale = TRUE)
      }
      jsonlite::write_json(autoscale_stats, path = file.path(limma_result_outdir, "autoscale_stats.json"))
    }
    if ("filter" %in% names(analysis)) {
      message("Filtering samples")
      for (criteria in analysis[["filter"]]) {
        samples_to_keep <- col_data[[criteria[["field"]]]] == criteria[["value"]]
        x <- x[samples_to_keep, ]
        col_data <- col_data[samples_to_keep,]
      }
    }
    design_formula <- as.formula(analysis[["formula"]])
    design <- model.matrix(
      design_formula,
      data = col_data %>% tibble::column_to_rownames("SampleID")
    )
    x_for_limma <- x[rownames(design), ]

    if (!all(rownames(x_for_limma) == rownames(design))) {
      stop("Error")
    }

    limma_fit <- limma::lmFit(t(x_for_limma), design = design, method= analysis[["fit_method"]])
    limma_fit <- limma::eBayes(limma_fit)
    saveRDS(limma_fit, file.path(limma_result_outdir, "limma_fit.rds"))
    de_results[[analysis_name]] <- list()
    for (result_contrast in analysis[["results"]]) {
      result_name <- result_contrast[["name"]]
      result_coef <- as.character(result_contrast[["coefficient"]])
      if (! result_coef %in% colnames(limma_fit)) {
        stop(
          paste0(
            "result ",
            result_coef,
            " not a coefficient. Coefs: ",
            paste0(colnames(limma_fit), collapse =", "))
        )
      }
      result_table <- limma::topTable(
        limma_fit,
        coef = result_coef,
        number = Inf,
        adjust.method = as.character(result_contrast[["padjust_method"]])
      )
      result_table <- tibble::rownames_to_column(result_table, "UniProtName")
      result_table <- tibble::as_tibble(result_table)
      result_table <- dplyr::select(
        result_table,
        "UniProtName" = "UniProtName",
        "baseMean" = "AveExpr",
        "log2FoldChange" = "logFC",
        "stat" = "t",
        "pvalue" = "P.Value",
        "padj" = "adj.P.Val"
      )
      if (!is.null(results_callback)) {
        result_table <- results_callback(result_table)
      }
      write.csv(
        x = result_table,
        file = file.path(limma_result_outdir, sprintf("%s.csv", result_name)),
        row.names = FALSE
      )
      saveRDS(
        result_table,
        file = file.path(limma_result_outdir, sprintf("%s.rds", result_name))
      )
      de_results[[analysis_name]][[result_name]] <- result_table
      ## Histogram of p-values
      png(filename = file.path(limma_result_outdir, sprintf("%s_hist_pval.png", result_name)))
      hist(result_table$pvalue[result_table$baseMean>1], xlab = "p-values (proteins with AveExpr > 1)", main = "Histogram")
      dev.off()
    }
  }
  de_results
}
