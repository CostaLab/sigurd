#'@import fastmatch
#'@param x vcf_path
`%fin%` <- function(x, table) {
  fmatch(x, table, nomatch = 0L) > 0L
}
