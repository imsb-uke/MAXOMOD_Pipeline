remove_ensembl_version <- function(ensembl_ids) {
  gsub(pattern = "(.*)\\.(.*)", replacement = "\\1", x = ensembl_ids)
}
