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
