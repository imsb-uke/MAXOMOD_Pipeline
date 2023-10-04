suppressPackageStartupMessages({
  requireNamespace("rlang")
  requireNamespace("ggplot2")
  requireNamespace("cowplot")
})

pc_scores_loadings_plot <- function(pca_model_info, pcx, pcy, color = NULL, color_values = NULL,
                                    loading_colname = "column", loading_label = "Feature") {
  slice_min_max <- function(...) {
    dplyr::bind_rows(
      dplyr::slice_min(...),
      dplyr::slice_max(...)
    )
  }

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
  ) +
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), 'cm')) 

  out_plot
}
