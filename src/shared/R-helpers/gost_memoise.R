source(here::here("src/shared/R-helpers/memoise_decorator.R"))

suppressPackageStartupMessages({
  requireNamespace("gprofiler2")
})


gost_memoise <- function(cache_dirname, ...) {

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


  fun <- memoise_decorator(
    cache_dirname = cache_dirname,
    cache_prefix = "gost",
    fun = .gost_robust
  )
  fun(...)
}
