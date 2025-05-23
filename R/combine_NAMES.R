#'combine_NAMES
#'@description
#'We combine two vectors of names.
#'@param x First vector of names.
#'@param y Second vector of names.
#'@export
combine_NAMES <- function(x, y) {
  shared_names <- intersect(x, y)
  combined_names <- c(x, setdiff(y, shared_names))
  return(combined_names)
}
