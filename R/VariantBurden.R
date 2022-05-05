#'@import Matrix SummarizedExperiment
#'@param se
VariantBurden <- function(se){
  burden <- colSums(assays(se)$fraction)
  colData(se)[,"Burden"] <- burden
  return(se)
}
