suppressPackageStartupMessages({
  requireNamespace("matrixStats")
})

get_top_n_variance <- function(matr, ntop = 500, tol = 1E-9) {
  rv <- matrixStats::colVars(matr)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # Filter low variance:
  select <- select[rv[select] > 1E-9]
  matr[, select]
}
