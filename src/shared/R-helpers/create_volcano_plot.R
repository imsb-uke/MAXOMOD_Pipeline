suppressPackageStartupMessages({
  requireNamespace("rlang")
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggrepel")
})

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
  } else {
    to_plot_known <- NULL
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
    if (!is.null(to_plot_known)) {
      to_plot_down_labels <- dplyr::anti_join(to_plot_down, to_plot_known, by = label_featurename)
    } else {
      to_plot_down_labels <- to_plot_down
    }
    gplt <- gplt + ggrepel::geom_text_repel(
      data = to_plot_down_labels,
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
    if (!is.null(to_plot_known)) {
      to_plot_up_labels <- dplyr::anti_join(to_plot_up, to_plot_known, by = label_featurename)
    } else {
      to_plot_up_labels <- to_plot_up
    }
    gplt <- gplt + ggrepel::geom_text_repel(
      data = to_plot_up_labels,
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
