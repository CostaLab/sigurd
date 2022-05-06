#'We combine two vectors of names.
#'@param x y
combine_NAMES <- function(x, y) {
  shared_names <- intersect(x, y)
  c(x, setdiff(y, shared_names))
}
