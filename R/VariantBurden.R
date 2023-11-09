#'VariantBurden
#'@description
#'Calculate the variant burden per cell.
#'We simply sum up the MAF values per cell.
#'@importFrom SummarizedExperiment assays colData
#'@importFrom Matrix colSums
#'@param se SummarizedExperiment object
#'@export
VariantBurden <- function(se){
  burden <- Matrix::colSums(SummarizedExperiment::assays(se)[["fraction"]])
  SummarizedExperiment::colData(se)[,"Burden"] <- burden
  return(se)
}
