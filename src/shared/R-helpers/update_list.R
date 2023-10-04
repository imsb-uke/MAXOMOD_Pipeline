update_list <- function(ours, to_update, recursiveness = 0L) {
  for (el in names(to_update)) {
    element <- to_update[[el]]
    if (recursiveness > 0L && !is.null(names(element))) {
      ours[[el]] <- update_list(ours[[el]], element, recursiveness = recursiveness - 1)
    } else {
      ours[[el]] <- element
    }
  }
  ours
}
