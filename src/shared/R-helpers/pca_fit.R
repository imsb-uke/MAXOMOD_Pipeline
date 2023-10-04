source(here::here("src/shared/R-helpers/get_top_n_variance.R"))

suppressPackageStartupMessages({
  requireNamespace("stats")
  requireNamespace("broom")
  requireNamespace("tidyr")
})


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
