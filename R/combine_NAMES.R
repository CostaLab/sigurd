#'@import BiocGenerics
#'@param se_somatic se_MT suffixes min_intersecting_cells
combine_NAMES <- function(x, y, ...) {
  shared_names <- intersect(x, y)
  c(x, setdiff(y, shared_names))
}
