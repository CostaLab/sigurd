#'VariantBurden
#'@description
#'Calculate the variant burden per cell.
#'We simply sum up the MAF values per cell.
#'@import Matrix SummarizedExperiment
#'@param se SummarizedExperiment object
#'@export
VariantBurden <- function(se){
  burden <- colSums(assays(se)[["fraction"]])
  colData(se)[,"Burden"] <- burden
  return(se)
}
