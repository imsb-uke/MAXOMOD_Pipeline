source(here::here("src/shared/R-helpers/memoise_decorator.R"))

suppressPackageStartupMessages({
  requireNamespace("clusterProfiler")
})


gseGO_memoise <- function(cache_dirname, ...) {

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

  fun <- memoise_decorator(
    cache_dirname = cache_dirname,
    cache_prefix = "gseGO",
    fun = .gseGO_robust
  )
  fun(...)
}
