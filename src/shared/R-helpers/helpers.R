suppressPackageStartupMessages({
  requireNamespace("argparse")
  requireNamespace("broom")
  requireNamespace("clusterProfiler")
  requireNamespace("cowplot")
  requireNamespace("digest")
  requireNamespace("dplyr")
  requireNamespace("gprofiler2")
  requireNamespace("matrixStats")
  requireNamespace("ggrepel")
  requireNamespace("rlang")
  requireNamespace("tidyr")
  requireNamespace("DESeq2")
#   library("readr")
#   library("tibble")
#   library("ggplot2")
#   library("preprocessCore")
   requireNamespace("yaml")
})

get_args <- function() {

  # create parser object
  parser <- argparse::ArgumentParser()

  parser$add_argument("--setting-key", type="character", default=NULL,
                      help="setting key from params.yaml to use (e.g. key1:subkey2")

  args <- parser$parse_args()
  args
}

get_settings <- function(...) {
  cmdline_args <- get_args()
  func_args <- list(...)

  final_args <- cmdline_args
  for (arg_name in names(func_args)) {
    final_args[[arg_name]] <- func_args[[arg_name]]
  }

  settings <- yaml::read_yaml("params.yaml")

  for (el in strsplit(final_args$setting_key, ":")[[1]]) {
    settings <- settings[[el]]
  }

  system_settings <- yaml::read_yaml("system_settings.yaml")

  for (el in strsplit(final_args$setting_key, ":")[[1]]) {
    system_settings <- system_settings[[el]]
  }

  for (el in names(system_settings)) {
    element <- system_settings[[el]]
    if (!is.null(names(element))) {
      for (ele in names(element)) {
        settings[[el]][[ele]] <- element[[ele]]
      }
    } else {
      settings[[el]] <- system_settings[[el]]
    }
  }
  settings
}


get_top_n_variance <- function(matr, ntop = 500, tol = 1E-9) {
  rv <- matrixStats::colVars(matr)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # Filter low variance:
  select <- select[rv[select] > 1E-9]
  matr[, select]
}

pca_fit <- function(matr, ntop = 500, cent = TRUE, scal = FALSE, scores_metadata = NULL) {
  pca <- stats::prcomp(get_top_n_variance(matr, ntop), center = cent, scale. = scal)
  var_percent <- pca$sdev^2 / sum(pca$sdev^2)
  pca_scores <- broom::augment(pca, matrix = "scores", data = scores_metadata)

  pca_loadings <- broom::tidy(pca, matrix = "loadings") %>%
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value")

  list(
    pca = pca,
    scores = pca_scores,
    loadings = pca_loadings,
    var_percent = var_percent
  )
}

head_tail <- function(x, n=6) {
  c(head(x, n), tail(x,n))
}

slice_min_max <- function(...) {
  dplyr::bind_rows(
    dplyr::slice_min(...),
    dplyr::slice_max(...)
  )
}


pc_scores_loadings_plot <- function(pca_model_info, pcx, pcy, color = NULL, color_values = NULL,
                                    loading_colname = "column", loading_label = "Feature") {
  xsym_scores <- rlang::sym(paste0(".fittedPC", as.character(pcx)))
  ysym_scores <- rlang::sym(paste0(".fittedPC", as.character(pcy)))
  loading_colname <- rlang::sym(loading_colname)
  gplt_scores <- ggplot2::ggplot(pca_model_info$scores)
  if (!is.null(color)) {
    colorsym <- rlang::sym(color)
    gplt_scores <- gplt_scores + ggplot2::geom_point(
      ggplot2::aes(x = !!xsym_scores, y = !!ysym_scores, color = !!colorsym, shape = Sex), size = 3
    )
  } else {
    gplt_scores <- gplt_scores + ggplot2::geom_point(
      ggplot2::aes(x = !!xsym_scores, y = !!ysym_scores, shape = Sex), size = 3
    )
  }
  if (!is.null(color_values)) {
    gplt_scores <- gplt_scores + ggplot2::scale_color_manual(values = color_values)
  }
  gplt_scores <- gplt_scores +
    ggplot2::labs(
      x = sprintf("PC%d (%2.0f%%)", pcx, 100*pca_model_info$var_percent[pcx]),
      y = sprintf("PC%d (%2.0f%%)", pcy, 100*pca_model_info$var_percent[pcy])) +
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")

  xsym_loads <- rlang::sym(sprintf("PC%d", pcx))
  gplt_loadx <- ggplot2::ggplot(slice_min_max(pca_model_info$loadings, !!xsym_loads, n = 10)) +
    ggplot2::geom_col(ggplot2::aes(x = forcats::fct_reorder(!!loading_colname, !!xsym_loads), y = !!xsym_loads)) +
    ggplot2::labs(x = loading_label, y = sprintf("PC%d loading", pcx)) +
    ggplot2::coord_flip()

  ysym_loads <- rlang::sym(sprintf("PC%d", pcy))
  gplt_loady <- ggplot2::ggplot(slice_min_max(pca_model_info$loadings, !!ysym_loads, n = 10)) +
    ggplot2::geom_col(ggplot2::aes(x = forcats::fct_reorder(!!loading_colname, !!ysym_loads), y = !!ysym_loads)) +
    ggplot2::labs(x = loading_label, y = sprintf("PC%d loading", pcy)) +
    ggplot2::coord_flip()


  out_plot <- cowplot::plot_grid(
    cowplot::plot_grid(NULL, gplt_scores, NULL, nrow = 3, ncol = 1,  rel_heights = c(1, 3, 1)),
    cowplot::plot_grid(gplt_loadx, gplt_loady, nrow = 2, ncol = 1),
    nrow = 1, ncol = 2
  )
  out_plot
}

cohort_to_factors <- function(cohort, cohort_factors) {
  out <- cohort
  for (cohort_factor in cohort_factors) {
    fac <- cohort_factor[["factor"]]
    if ("levels" %in% names(cohort_factor)) {
      out[[fac]] <- factor(out[[fac]], levels = as.character(cohort_factor[["levels"]]))
    } else {
      out[[fac]] <- factor(out[[fac]])
    }
  }
  out
}

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
          mean= mean(col_data[[scale_var]], na.rm = TRUE),
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
      ## Histogram of p-values
      png(filename = file.path(deseq2_result_outdir, sprintf("%s_hist_pval.png", result_contrast[["name"]])))
      hist(r$pvalue[r$baseMean>1], xlab = "p-values (genes with baseMean > 1)", main = "Histogram")
      dev.off()
    }
  }
}


create_volcano_plot <- function(
  df,
  padj_thres = 0.05,
  up_log2fc_thres = log2(1.5),
  down_log2fc_thres = -log2(1.5),
  padj_thres_loose=0.1,
  cap_minuslog10pval = 10,
  label_featurename = NULL,
  label_featurefiltername = NULL,
  label_featurefiltervalues=NULL) {
  # force so ggplot works (avoid lazy evaluation issues)
  force(padj_thres)
  force(up_log2fc_thres)
  force(down_log2fc_thres)
  force(padj_thres_loose)
  force(label_featurename)
  force(cap_minuslog10pval)
  force(label_featurefiltername)
  force(label_featurefiltervalues)
  if (!is.null(label_featurename)) {
    label_featurenamesym <- rlang::sym(label_featurename)
  }
  if (!is.null(label_featurefiltername)) {
    label_featurefilternamesym <- rlang::sym(label_featurefiltername)
  }
  to_plot <- df
  to_plot$capped_minuslog10pval <- pmin(-log10(to_plot$pvalue), cap_minuslog10pval)
  to_plot_down <- dplyr::filter(to_plot, padj < !!padj_thres, log2FoldChange < !!down_log2fc_thres)
  to_plot_up <- dplyr::filter(to_plot, padj < !!padj_thres, log2FoldChange > !!up_log2fc_thres)
  to_plot_down_loose <- dplyr::filter(to_plot, padj < !!padj_thres_loose, padj > !!padj_thres, log2FoldChange < !!down_log2fc_thres)
  to_plot_up_loose <- dplyr::filter(to_plot, padj < !!padj_thres_loose, padj > !!padj_thres, log2FoldChange > !!up_log2fc_thres)
  if (!is.null(label_featurefiltername) && !is.null(label_featurefiltervalues)) {
    to_plot_known <- dplyr::filter(to_plot, !!label_featurefilternamesym %in% label_featurefiltervalues)
  }
  gplt <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = to_plot,
      mapping = ggplot2::aes(
        x = log2FoldChange,
        y = capped_minuslog10pval
      ),
      size=0.4, color = "gray") +
    ggplot2::geom_hline(yintercept = -log10(padj_thres), color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = down_log2fc_thres, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = up_log2fc_thres, color = "red", linetype = "dashed") +
    ggplot2::geom_point(
      data = to_plot_down,
      mapping = ggplot2::aes(x = log2FoldChange, y = capped_minuslog10pval),
      size=1,
      color = "blue"
    )
    if (!is.null(label_featurename)) {
      gplt <- gplt + ggrepel::geom_text_repel(
        data = to_plot_down,
        mapping = ggplot2::aes(x = log2FoldChange, y = capped_minuslog10pval, label = !!label_featurenamesym)
      )
    }
    gplt <- gplt + ggplot2::geom_point(
      data = to_plot_up,
      mapping = ggplot2::aes(x = log2FoldChange, y = capped_minuslog10pval),
      size=1,
      color = "red"
      )
    if (!is.null(label_featurename)) {
      gplt <- gplt + ggrepel::geom_text_repel(
        data = to_plot_up,
        mapping = ggplot2::aes(x = log2FoldChange, y = capped_minuslog10pval, label = !!label_featurenamesym)
      )
    }
    gplt <- gplt + ggplot2::geom_point(
      data = to_plot_down_loose,
      mapping = ggplot2::aes(x = log2FoldChange, y = capped_minuslog10pval),
      size=1, color = "#9696ff"
    ) +
    ggplot2::geom_point(
      data = to_plot_up_loose,
      mapping = ggplot2::aes(x = log2FoldChange, y = capped_minuslog10pval),
      size=1,
      color = "#ff9696"
    )
    if (!is.null(label_featurename) &&
        !is.null(label_featurefiltername) &&
        !is.null(label_featurefiltervalues)
      ) {
      gplt <- gplt + ggrepel::geom_text_repel(
        data = to_plot_known,
        mapping = ggplot2::aes(x = log2FoldChange, y = capped_minuslog10pval, label = !!label_featurenamesym)
      )
    }
    gplt <- gplt +
      ggplot2::labs(x = "log2FoldChange", y = "-log10(pvalue)") +
      ggplot2::theme_bw()
    gplt
}

memoise_decorator <- function(cache_dirname, cache_prefix, fun) {
  cache_dirname <- force(cache_dirname)
  cache_prefix <- force(cache_prefix)
  dir.create(cache_dirname, showWarnings = FALSE, recursive = TRUE)
  function(...) {
    args <- list(...)
    hash_value <- digest::digest(args)
    cache_fn <- file.path(cache_dirname, paste0(cache_prefix, "_", hash_value, ".rds"))
    if (file.exists(cache_fn)) {
      out <- readRDS(cache_fn)
    } else {
      out <- fun(...)
      saveRDS(out, file = cache_fn)
    }
    out
  }
}

.gost_robust <- function(...) {
  args <- list(...)
  q <- args$query
  if (length(q) > 0) {
    Sys.sleep(1)
    out <- gprofiler2::gost(...)
  } else {
    out <- NULL
  }
  if (is.null(out)) {
    out <- list(
      result = structure(
        list(
          query = character(0),
          significant = logical(0),
          p_value = numeric(0),
          term_size = integer(0),
          query_size = integer(0),
          intersection_size = integer(0),
          precision = numeric(0),
          recall = numeric(0),
          term_id = character(0),
          source = character(0),
          term_name = character(0),
          effective_domain_size = integer(0),
          source_order = integer(0),
          parents = list(),
          evidence_codes = character(0),
          intersection = character(0)
        ),
        class = "data.frame",
        row.names = integer(0)
      ),
      meta = NULL
    )
  }
  out
}

gost_memoise <- function(cache_dirname, ...) {
  fun <- memoise_decorator(
    cache_dirname = cache_dirname,
    cache_prefix = "gost",
    fun = .gost_robust
  )
  fun(...)
}

.enrichGO_robust <- function(...) {
  args <- list(...)
  gene <- args$gene
  if (length(gene) == 0) {
    out <- NULL
  } else {
    Sys.sleep(1)
    out <- clusterProfiler::enrichGO(...)
  }
  out
}

enrichGO_memoise <- function(cache_dirname, ...) {
  fun <- memoise_decorator(
    cache_dirname = cache_dirname,
    cache_prefix = "enrichGO",
    fun = .enrichGO_robust
  )
  fun(...)
}


.gseGO_robust <- function(...) {
  args <- list(...)
  geneList <- args$geneList
  if (length(geneList) == 0) {
    out <- NULL
  } else {
    Sys.sleep(1)
    out <- clusterProfiler::gseGO(...)
  }
  out
}


gseGO_memoise <- function(cache_dirname, ...) {
  fun <- memoise_decorator(
    cache_dirname = cache_dirname,
    cache_prefix = "gseGO",
    fun = .gseGO_robust
  )
  fun(...)
}


.gseKEGG_robust <- function(...) {
  args <- list(...)
  geneList <- args$geneList
  if (length(geneList) == 0) {
    out <- NULL
  } else {
    Sys.sleep(1)
    out <- clusterProfiler::gseKEGG(...)
  }
  out
}


gseKEGG_memoise <- function(cache_dirname, ...) {
  fun <- memoise_decorator(
    cache_dirname = cache_dirname,
    cache_prefix = "gseKEGG",
    fun = .gseKEGG_robust
  )
  fun(...)
}

remove_ensembl_version <- function(ensembl_ids) {
  gsub(pattern = "(.*)\\.(.*)", replacement = "\\1", x = ensembl_ids)
}
