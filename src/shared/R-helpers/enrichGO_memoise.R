source(here::here("src/shared/R-helpers/memoise_decorator.R"))

suppressPackageStartupMessages({
  requireNamespace("clusterProfiler")
})


enrichGO_memoise <- function(cache_dirname, ...) {
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

  fun <- memoise_decorator(
    cache_dirname = cache_dirname,
    cache_prefix = "enrichGO",
    fun = .enrichGO_robust
  )
  fun(...)
}
