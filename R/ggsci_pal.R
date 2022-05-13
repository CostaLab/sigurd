#'Function to return colours from a ggsci palette.
#'@import assertthat ggsci glue
#'@param option Your colour palette of choice.
#'@details
#'The function returns a colour palette from ggsci.
#'Options are:
#'aaas: 10 d3: 10 futurama: 12 gsea: 12
#'igv: 51 jama: 7 jco: 10 npg: 10 lancet: 9 
#'locuszoom: 7 material: 10 nejm: 8 
#'rickandmorty: 12 simpsons: 16 
#'startrek: 7 tron: 7 uchicago: 9 
#'ucscgb: 26
#'@export
ggsci_pal <- function(option, ...){
  func_name = glue("pal_{option}")
  func_call = glue('{func_name}(...)')
  assertthat::assert_that(func_name %in% ls("package:ggsci"))
  return(eval(parse(text=func_call)))
}
