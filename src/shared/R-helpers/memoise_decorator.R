suppressPackageStartupMessages({
  requireNamespace("digest")
})

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
